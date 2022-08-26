require(tidyverse)
require(scales)

a = read_tsv("v52.2_HO_public.anno") %>%
    mutate(lat=as.numeric(Lat.), long=as.numeric(Long.)) %>%
    rename(year=6, age=9) %>%
    mutate(year=as.numeric(year), age=as.numeric(age)) %>%
    filter(age>0, year<=2022, !is.na(year), !is.na(lat), !is.na(long)) %>%
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

make_fig <- function(a, fname){
P2 = ggplot() + 
    geom_sf(data=world, color=NA) + 
    theme_minimal() + 
    geom_point(data=a, aes(color=age, x=long, y=lat)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    #theme(legend.position = 'none') + 
    xlab(NULL) + ylab(NULL) +
    coord_sf(xlim=c(-170, 170), ylim=c(-55, 70)) +
    scale_color_viridis_c(breaks=c(0000, 10000, 20000),
    limits=c(0000,20000), oob=squish)

ggsave(fname, P2, width=7, height=3)
}

a %>% filter(year <= 2010) %>% make_fig(fname='figures/ancient_dna_till2010.png')
a %>% filter(year <= 2014) %>% make_fig(fname='figures/ancient_dna_till2014.png')
a %>% filter(year <= 2016) %>% make_fig(fname='figures/ancient_dna_till2016.png')
a %>% filter(year <= 2022) %>% make_fig(fname='figures/ancient_dna_till2022.png')

a %>% filter(age > 40000) %>% make_fig(fname='figures/ancient_dna_older_than40k.png')


a %>% ggplot(aes(x=age)) + geom_histogram(binwidth=1000) -> P3
ggsave('figures/ancient_dna_by_year.png', P3, width=5, height=3)
