library("rnaturalearth")
library("rnaturalearthdata")
library(ggrepel)
world <- ne_countries(scale = "medium", returnclass = "sf")
source("fscripts.R")
library(ggpubr)
BTHEME = theme_classic() + theme(legend.position="none")


anno = read_tsv("v52.2_HO_public.anno") %>%
    select(iid=2,
           age=6,
           pop=8,
           lat=11,
           long=12,
           src=13,
           cov=14) %>%
    mutate(lat = as.numeric(lat),
           long=as.numeric(long),
           age = as.numeric(age),
           cov = as.numeric(cov))

ind_eu = read_table("subdata/westeurasian1.ind", col_names=c("iid", "sex", "pop"))
ind_world = read_table("subdata/worldfoci2.ind", col_names=c("iid", "sex", "pop"))

pops_eu = ind_eu %>% left_join(anno) %>% select(pop, lat, long) %>% 
    distinct() %>%
    group_by(pop) %>% summarize_all(mean)
pops_world = ind_world %>% left_join(anno) %>% select(pop, lat, long) %>% distinct() %>%
    group_by(pop) %>% summarize_all(mean)


#mapp

plot_map_eu <- function(){
    P1 = ggplot() +
        geom_sf(data=world, color=NA) +
        theme_minimal() +
        geom_point(data=pops_eu, aes(x=long, y=lat)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(legend.position = 'none') +
        xlab(NULL) + ylab(NULL) +
        coord_sf(xlim=range(pops_eu$long,na.rm=T), ylim=range(pops_eu$lat, na.rm=T)) +
        scale_color_viridis_c(trans='log')
}

P1 = plot_map_eu() + 
    geom_text_repel(data=pops_eu, aes(x=long, y=lat, label=pop), 
                    size=2, segment.color='darkgrey') 
ggsave("figures/sample_map_eu.png", P1, width=4, height=3)
ggsave("figures/sample_map_eu2.png", plot_map_eu(), width=4, height=3)

plot_map_world <- function(){
    P2 = ggplot() +
        geom_sf(data=world, color=NA) +
        theme_minimal() +
        geom_point(data=pops_world, aes(x=long, y=lat)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(legend.position = 'none') +
        xlab(NULL) + ylab(NULL) +
        coord_sf(xlim=range(pops_world$long,na.rm=T), ylim=range(pops_world$lat, na.rm=T)) +
        scale_color_viridis_c(trans='log')
}

P2 = plot_map_world() + 
    geom_text_repel(data=pops_world, aes(x=long, y=lat, label=pop), size=2,
                    segment.color='darkgrey') 
ggsave("figures/sample_map_world.png", P2, width=7, height=3)
ggsave("figures/sample_map_world2.png", plot_map_world(), width=4, height=3)


