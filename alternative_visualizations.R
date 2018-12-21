my_frame_ultra = my_frame
my_frame_ultra%>%filter(cromo_i == cromo_j)
my_frame_ultra%>%filter(cromo_i != cromo_j, cromo_i==18)

my_frame_ultra = 
my_frame_ultra%>%
  mutate(ultra.intra.inter = ifelse(test = cromo_i!=cromo_j, 
                                    yes = "inter", 
                                    no = ifelse(
                                      test = cromo_i%in%c(1,2,3,19,20,"X"), 
                                      yes = "metacentric_intra",
                                      no = ifelse(test = cromo_i%in%c(13,14,15,21,22), 
                                                  yes = "acrocentric_intra",
                                                  no = "submetacentric_intra"
                                                    )
                                      )
                                    )
         )
my_frame_ultra%>%head

p = ggplot(data = my_frame_ultra, mapping = aes(signFactor, 
                                          log(changeFactor), 
                                          colour = factor(ultra.intra.inter),
                                          label  = paste0(cromo_i, ";", cromo_j)
)
)
p = p + geom_point()
p = p + geom_text(check_overlap=FALSE, size = 3,show.legend = FALSE)
p = p + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p = p + labs(color='-chromosomal interaction') 
plot(p)
# ggsave(filename = "results/scatter_signchangefactors_basalVhealth2.pdf", 
#        height = 8.5, 
#        width = 14, 
#        units = "in")


my_frame_ultra%>%filter(ultra.intra.inter == "inter")


my_frame_ultra%>%filter(!(ultra.intra.inter == "inter" & cromo_i==18))%>%
ggplot(mapping = aes(signFactor, 
                                                log(changeFactor), 
                                                colour = factor(ultra.intra.inter),
                                                label  = paste0(cromo_i, ";", cromo_j)
)
)+ geom_point()+ geom_text(check_overlap=FALSE, size = 3,show.legend = FALSE) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+ labs(color='-chromosomal interaction') 
 ggsave(filename = "results/scatter_centricinfo.pdf", 
        height = 8.5, 
        width = 14, 
        units = "in")
 
 my_frame_ex = 
 my_frame%>%mutate(cromo_i  = factor(cromo_i, levels = c(1:22, "X"), ordered = TRUE),
                   cromo_j  = factor(cromo_j, levels = c(1:22, "X"), ordered = TRUE))
 
 my_frame_ex =
 my_frame_ex%>%filter(!(intra.inter == "inter" & cromo_i > cromo_j))
 
 
 
 