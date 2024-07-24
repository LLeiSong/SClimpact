#' @title Hexify the continents.
#' @description Take the global continents map and make hex grids.
#' @param res (`numeric`) The resolution in meter.
#' @param crs (`crs` or `numeric`) The projected coordinate system. 
#' It can be an object of `crs` or EPSG number. Only projected coordinates 
#' systems are allowed to make the area the same for all hex grids.
#' @param flat_topped (`logical`) A parameter for \code{\link{st_make_grid}}.
#' If TRUE generate flat topped hexagons, else generate pointy topped.
#' @return The hex grids with columns of country name, ISOs and continent name.
#'
#' @seealso
#' \code{\link{st_make_grid}}, \code{\link{ne_countries}}
#'
#' @details
#' Convert the terrestrial continents into hex grids. Note that the map will be
#' shifted to keep Russia together. So it will be problematic if the user 
#' directly transform to geographic coordinate system.
#'
#' @importFrom sf st_make_valid st_transform st_union st_make_grid st_as_sf st_join
#' @importFrom dplyr filter select rename slice mutate
#' @importFrom rnaturalearth ne_countries
#' @importFrom lwgeom st_wrap_x
#' @export
#' @examples
#' # Get the hex grids with 500km at pseudo Mercator projection system
#' hex_grids <- hexify_continents(res = 100000, crs = 3857)

hexify_continents <- function(res = 5000,
                              crs = 3857,
                              flat_topped = FALSE){
    # Check sf setting
    org_val <- sf_use_s2()
    sf_use_s2(FALSE)
    
    # Query the boundary with small scale
    terr_bry <- ne_countries(scale = 110) %>% 
        filter(continent != "Antarctica") %>% 
        filter(continent != "Seven seas (open ocean)") %>% 
        st_transform(crs) %>% select(formal_en, iso_a2, iso_a3, continent)
    
    bbox_orig <- st_bbox(terr_bry)
    
    # Move the min x from -20037508 to -18757508 to put Russia together
    terr_bry <- terr_bry %>% filter(iso_a2 != "RU") %>% 
        rbind(st_wrap_x(terr_bry %>% filter(iso_a2 == "RU"), 
                        -18757508, bbox_orig[3] - bbox_orig[1]))
        
    # Mosaic Russia
    ru <- terr_bry %>% filter(iso_a2 == "RU")
    ru <- ru %>% mutate(geometry = st_union(st_buffer(ru, 100)) %>% st_buffer(-100))
    terr_bry <- terr_bry %>% filter(iso_a2 != "RU") %>% 
        rbind(ru)
    
    # Union up the boundary to make grids
    bry <- terr_bry %>% st_make_valid() %>% st_union()
    
    # Index the disconnect lands
    lands <- st_cast(bry, "POLYGON") %>% st_as_sf() %>% 
        mutate(land_id = 1:nrow(.)) %>% group_by(land_id) %>% 
        summarise(geometry = st_union(x))
    
    # Merge Russian as the whole
    lands <- lands %>% mutate(land_id = ifelse(land_id == 18, 40, land_id)) %>% 
        mutate(land_id = ifelse(land_id == 21, 115, land_id)) %>% 
        group_by(land_id) %>% summarise(geometry = st_union(geometry))
    
    # Remove the island that is smaller than the grid size
    lands <- lands %>% mutate(area = st_area(.)) %>% 
        filter(area >= units::set_units(sqrt(3) * res^2 / 2, "m2")) %>% 
        select(-area) %>% mutate(land_id = 1:nrow(.))
    
    # Make hex grids
    hex_grids <- st_make_grid(bry, res, crs = crs, 
                              square = FALSE, flat_topped = flat_topped)
    
    if (length(hex_grids) > 160000) {
        # Chunk up to avoid RAM crash
        hex_grids <- hex_grids %>% st_as_sf() %>% rename(geometry = x) %>% 
            slice(unique(unlist(st_intersects(terr_bry, .)))) %>% 
            mutate(gid = 1:nrow(.))
        
        # Split into groups
        chunklength <- 1e4
        grps <- split(1:nrow(hex_grids), 
                      ceiling(seq_along(1:nrow(hex_grids)) / chunklength))
        hex_grids <- do.call(rbind, pblapply(grps, function(ids){
            # message(ids)
            hex_grids %>% slice(ids) %>% 
                # Add a tiny buffer to avoid invalid intersects
                st_join(st_buffer(terr_bry, -res/100000), largest = TRUE) %>% 
                st_join(st_buffer(lands, -res/100000), largest = TRUE)
        })) %>% arrange(gid)
    } else {
        # Direct cal for light data
        hex_grids <- hex_grids %>% st_as_sf() %>% rename(geometry = x) %>% 
            slice(unique(unlist(st_intersects(terr_bry, .)))) %>% 
            mutate(gid = 1:nrow(.)) %>% 
            # Add a tiny buffer to avoid invalid intersects
            st_join(st_buffer(terr_bry, -res/100000), largest = TRUE) %>% 
            st_join(st_buffer(lands, -res/100000), largest = TRUE)
    }
    
    # Clean the grids and re-index
    hex_grids <- hex_grids %>% filter(!is.na(land_id)) %>% mutate(gid = 1:nrow(.))
    
    # # Shift the bounding back
    # ## Define the line for the cut
    # ls <- rbind(c(bbox_orig[3], st_bbox(hex_grids)[2]), 
    #             c(bbox_orig[3], st_bbox(hex_grids)[4]))
    # ls <- st_linestring(ls)
    # ls <- st_sfc(ls, crs = crs)
    # ids <- unique(unlist(st_intersects(ls, hex_grids)))
    # grids_ls <- hex_grids %>% slice(ids)
    # 
    # ## Rearrange the grids
    # hex_grids <- st_wrap_x(
    #     hex_grids[-ids, ], bbox_orig[3], 
    #     -(st_bbox(hex_grids)[3] - st_bbox(hex_grids)[1]))
    # 
    # ## Remove partial grids
    # hex_grids <- rbind(hex_grids, grids_ls)
    
    # Change sf setting back
    sf_use_s2(org_val)
    
    # Return the object
    return(hex_grids)
}

# # Generate useful hex grids
# for (rs in c(200, 100, 50, 10, 5)){
#     message(sprintf("Generate hex grids with res: %s", rs))
#     hex_grids <- hexify_continents(res = rs * 1000)
#     st_write(hex_grids, sprintf("data/geoms/hex_grids_%skm_global.geojson", rs))
# }
# # Post clean up
# rm(list=ls()); gc()
