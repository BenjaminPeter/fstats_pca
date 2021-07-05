require(tidyverse)
a  = read_tsv("v44.3_1240K_public.anno") 

z = a %>% rename(year=6,age=10,coverage=20) %>% 
    mutate(year=as.numeric(year), age=as.numeric(age)) %>% 
    mutate(lat=as.numeric(Lat.), long=as.numeric(Long.)) %>%
    filter(age>0)

P = z %>% filter(!is.na(year), age>0) %>% group_by(year) %>% 
    tally %>% mutate(n=cumsum(n)) %>% 
    ungroup %>%
    ggplot(aes(x=as.factor(year), y=n)) + geom_col() + 
    xlab("Publication Year") + ylab("# samples")
ggsave("figures/ancient_dna_by_year.png", P, width=5, height=3)


P2= z %>% mutate(coverage=as.numeric(coverage), coverage=pmin(coverage,10)) %>% 
    ggplot(aes(x=as.numeric(coverage))) + geom_histogram(bins=100) 
#    xlab=("Coverage") + ylab=("# samples")
ggsave("figures/ancient_dna_coverage.png", P2, width=5, height=3)


require(rnaturalearth)
require(rnaturalearthdata)
world <- ne_countries(scale = "medium", returnclass = "sf")
P3 = ggplot() + 
    geom_sf(data=world, color=NA) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point(aes(color=age, x=long, y=lat), data=z) + 
    scale_color_viridis_c(trans='log') +
    xlab("") + ylab("") + 
    theme(legend.pos='none') + 
    coord_sf(xlim=c(-170, 160), ylim=c(-50,75))
ggsave("figures/ancient_dna_map.png", P3, width=7, height=4, dpi=150)
