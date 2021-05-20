require(tidyverse)

a = read_tsv("v44.3_1240K_public.anno") %>%
    mutate(lat=as.numeric(Lat.), long=as.numeric(Long.)) %>%
    rename(year=6, age=10) %>%
    mutate(year=as.numeric(year), age=as.numeric(age)) %>%
    filter(age>0, !is.na(year), !is.na(lat), !is.na(long)) %>%
    select(lat, long, year, age)

df_year = a %>% 
    mutate(year=as.factor(year)) %>% 
    group_by(year) %>% 
    tally %>% 
    mutate(n=cumsum(n))

P = df_year %>% ggplot(aes(x=year, y=n)) + geom_col()
ggsave("figures/ancient_dna_by_year.png", P, width=5, height=3)

library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")

P2 = ggplot() + 
    geom_sf(data=world, color=NA) + 
    theme_minimal() + 
    geom_point(data=a, aes(color=age, x=long, y=lat)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position = 'none') + 
    xlab(NULL) + ylab(NULL) +
    coord_sf(xlim=c(-170, 170), ylim=c(-55, 70)) +
    scale_color_viridis_c(trans='log')

ggsave("figures/ancient_dna_map.png", P2, width=7, height=3)
