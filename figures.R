# Original manuscript layout, using the fitted objects from replicate.R.
library(ggplot2)
library(ggrepel)
library(patchwork)
library(sf)
type_colours <- c(RPPN="orange","Urban park"="orangered2")
export_plot <- function(number, plot, width, height) {
  for(ext in c("png","pdf"))
    ggsave(sprintf("results/Figure_%d.%s",number,ext),plot=plot,
           width=width,height=height,units="in",dpi=300,bg="white")
}
# Cache boundary geometry only. All fitted landscape models use the supplied
# 600 m land-cover table and need no external raster or geographic service.
boundary_file <- "data/curitiba_boundaries.rds"
if(file.exists(boundary_file)) {
  bounds <- readRDS(boundary_file)
} else {
  bairros <- geobr::read_neighborhood(year=2010,simplified=TRUE)
  bounds <- list(neighbourhoods=st_transform(bairros[bairros$code_muni==4106902,],31982),
                 city=st_transform(geobr::read_municipality(code_muni=4106902,year=2020),31982))
  stopifnot(nrow(bounds$neighbourhoods)>0,nrow(bounds$city)==1)
  saveRDS(bounds,boundary_file)
}
pts <- st_as_sf(transform(sites,Site=rownames(sites)),coords=c("easting_m","northing_m"),crs=31982)
centre <- st_sfc(st_point(c(673648.4910977068,7186491.014442006)),crs=31982)
location <- cbind(st_drop_geometry(pts),st_coordinates(pts))
map <- ggplot() +
  geom_sf(data=bounds$neighbourhoods,fill="grey94",colour="grey70",linewidth=.18) +
  geom_sf(data=bounds$city,fill=NA,colour="black",linewidth=.50) +
  geom_sf(data=st_buffer(pts,600),fill=NA,colour="grey40",linewidth=.30) +
  geom_sf(data=pts,aes(colour=Type),size=2.4) +
  geom_sf(data=centre,shape=3,size=2.6) +
  geom_text_repel(data=location,aes(X,Y,label=Site),size=2.7,seed=12,max.overlaps=Inf) +
  scale_colour_manual(values=type_colours,name="Site type") +
  coord_sf(crs=st_crs(31982),datum=st_crs(4674)) +
  theme_bw() + theme(legend.position="bottom")
export_plot(1,map,7,7)

scr <- data.frame(fit$x[,1:2],Site=site)
hull <- do.call(rbind,lapply(split(scr,scr$Site),function(g) g[chull(g$PC1,g$PC2),]))
vec <- data.frame(fit$rotation[,1:2]*8,element=rownames(fit$rotation))
axis_name <- function(j) sprintf("PC%d (%.1f%%)",j,100*fit$sdev[j]^2/sum(fit$sdev^2))
pca_panel <- ggplot(scr,aes(PC1,PC2,colour=Site)) +
 geom_polygon(data=hull,aes(fill=Site,group=Site),alpha=.10,colour=NA,show.legend=FALSE) +
 geom_point(size=2.4,alpha=.88) +
 geom_segment(data=vec,aes(x=0,y=0,xend=PC1,yend=PC2),inherit.aes=FALSE,
              linewidth=.3,arrow=grid::arrow(length=grid::unit(.11,"cm"))) +
 geom_text_repel(data=vec,aes(PC1,PC2,label=element),inherit.aes=FALSE,
                 size=2.8,seed=12,min.segment.length=Inf) +
 scale_colour_viridis_d(option="H") + scale_fill_viridis_d(option="H") +
 labs(x=axis_name(1),y=axis_name(2)) + theme_bw() +
 theme(panel.grid=element_blank(),legend.key=element_blank(),legend.text=element_text(size=7))
# Display order and original reporting units follow the submitted heatmap.
keys <- c("P","K","Mn","Al","Fe","Zn","Co","Cr","Cu")
unit_labels <- setNames(names(raw),colnames(chem))
heat <- rbind(ct,transform(ct,a=b,b=a))
heat$a <- factor(heat$a,levels=keys);heat$b <- factor(heat$b,levels=keys)
heat$star <- ifelse(heat$q<.001,"***",ifelse(heat$q<.01,"**",ifelse(heat$q<.05,"*","")))
hm <- ggplot(heat,aes(a,b,fill=rho)) + geom_tile(colour="white",linewidth=.15) +
 geom_text(aes(label=star),size=3) +
 scale_fill_gradient2(low="royalblue3",mid="white",high="firebrick3",
                      limits=c(-1,1),name="Spearman\nrho") +
 scale_x_discrete(labels=unit_labels)+scale_y_discrete(labels=unit_labels)+
 coord_equal()+theme_classic()+labs(x=NULL,y=NULL)+
 theme(axis.text=element_text(size=7),axis.text.x=element_text(angle=45,hjust=1))
export_plot(2,wrap_plots(pca_panel,hm,nrow=1,widths=c(1.12,1)),10,5.54)

sd_plot <- data.frame(sites,sc,Site=rownames(sites))
dist_panel <- function(pc) {
 ggplot(sd_plot,aes(distance_km,.data[[pc]])) +
 geom_smooth(method="lm",formula=y~x,se=TRUE,colour="black",fill="grey80",
             alpha=.55,linewidth=.7) +
 geom_point(aes(colour=Type),size=3.3)+
 geom_text_repel(aes(label=Site),size=2.8,seed=12,max.overlaps=Inf)+
 scale_colour_manual(values=type_colours,name="Site type")+
 labs(x="Distance to Curitiba city centre (km)",y=paste(pc,"(site mean)"))+
 theme_bw()+theme(panel.grid=element_blank(),legend.key=element_blank())
}
export_plot(3,wrap_plots(dist_panel("PC1"),dist_panel("PC2"),nrow=1,guides="collect"),10,5)
