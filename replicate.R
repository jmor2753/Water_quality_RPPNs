# Run from the extracted folder: Rscript replicate.R
library(vegan)
dir.create("results",showWarnings=FALSE)
out <- function(x,name) write.csv(x,file.path("results",paste0(name,".csv")))
d <- read.csv("data/chemistry.csv",check.names=FALSE,stringsAsFactors=FALSE)
sites <- read.csv("data/sites.csv",row.names=1,check.names=FALSE)
site <- factor(d$site)
stopifnot(nrow(d)==36,all(table(site)==3),setequal(levels(site),rownames(sites)))
raw <- d[, -(1:2)]
flags <- sapply(raw,function(x) x %in% c("<LOQ","<LQO"))
chem <- sapply(raw,function(x) as.numeric(replace(x,x %in% c("<LOQ","<LQO"), NA)))
stopifnot(all(is.na(chem)==flags))
mg <- grepl("(mg/L)",names(raw),fixed=TRUE)
chem[,mg] <- chem[,mg]*1000
colnames(chem) <- colnames(flags) <- sub(" .*","",names(raw))
rownames(chem) <- d$sample_id
out(colSums(flags),"censor_counts")
chem[which(chem<0)] <- 0
chem[flags] <- .001
el <- c("Al","Cr","Cu","Fe","K","Mn","P","Zn")
base <- chem[,el]
out(chem,"analysis_concentrations_ug_L")
perm <- function(n) as.matrix(read.csv(gzfile(sprintf("data/permutations_%d.csv.gz",n)),header=FALSE))
P36 <- perm(36); P12 <- perm(12)
ord <- function(x) {
  fit <- prcomp(log10(1+x),scale.=TRUE)
  direction <- ifelse(fit$rotation["Fe",]<0,-1,1)
  direction[2] <- -direction[2] # Match the submitted PC2 orientation (Zn positive).
  fit$x <- sweep(fit$x,2,direction,"*")
  fit$rotation <- sweep(fit$rotation,2,direction,"*")
  fit
}
fit <- ord(base)
out(fit$x,"scores");out(fit$rotation,"loadings")
out(fit$sdev^2/sum(fit$sdev^2),"variance")
# Permutations reassign complete sample vectors to fixed groups; preserve the
# original archive's alphabetical site order and within-site input order.
idx <- order(site,seq_along(site))
scenario <- function(x) {
  y <- scale(log10(1+x))[idx,,drop=FALSE]
  group <- site[idx]
  parts <- adonis2(y~group,method="euclidean",permutations=P36)
  f <- ord(x)
  a <- fit$x[,1:2]; b <- f$x[,1:2]
  r <- sum(svd(crossprod(a,b),nu=0,nv=0)$d)/sqrt(sum(a^2)*sum(b^2))
  c(PC1=unname(f$sdev[1]^2/sum(f$sdev^2)),PC2=unname(f$sdev[2]^2/sum(f$sdev^2)),Procrustes_r=r,site_R2=parts$R2[1],F=parts$F[1],p=parts[1,"Pr(>F)"])
}
sens <- list(primary=scenario(base))
for (omit in list("Cr","Zn",c("Cr","Zn"))) sens[[paste("omit",paste(omit,collapse="+"))]] <- scenario(base[,setdiff(el,omit)])
for(frac in c(0,.25,.5,1)) {
  x <- base
  for(k in c("Cr","Zn")) x[flags[,k],k] <- frac*min(base[!flags[,k],k])
  sens[[paste("minimum times",frac)]] <- scenario(x)
}
out(do.call(rbind,sens),"sensitivity")
z <- scale(log10(1+base))
out(as.data.frame(anova(betadisper(dist(z),site,type="median"))),"dispersion")
means <- function(x) rowsum(x,site)/3
out(means(base),"Table_1_site_means_ug_L")
active <- chem[,apply(chem,2,sd)>0]
corr <- function(x) {
  pairs <- combn(colnames(x),2,simplify=FALSE)
  tab <- do.call(rbind,lapply(pairs,function(ab) {
    test <- cor.test(x[,ab[1]],x[,ab[2]],method="spearman",exact=FALSE)
    data.frame(a=ab[1],b=ab[2],rho=unname(test$estimate),p=test$p.value)
  }))
  tab$q <- p.adjust(tab$p,"BH");tab
}
ct <- corr(active);out(ct,"Table_S2_spearman_samples");out(corr(means(active)),"site_correlations")
sites$distance_km <- sqrt((sites$easting_m-673648.4910977068)^2+(sites$northing_m-7186491.014442006)^2)/1000
sc <- means(fit$x)[rownames(sites),1:2]
response <- means(z)[rownames(sites),]
reg <- list(); rd <- list()
for(k in c("distance_km","urban_percent","forest_percent")) {
  x <- sites[[k]]
  for(j in 1:2) {
    f <- summary(lm(sc[,j]~x))
    reg[[length(reg)+1]] <- data.frame(predictor=k,axis=colnames(sc)[j],slope=coef(f)[2,1],intercept=coef(f)[1,1],R2=f$r.squared,p=coef(f)[2,4])
  }
  mod <- rda(response~x)
  # Permuting predictors against fixed responses equals permuting responses
  # by the inverse permutations. Retain the archived randomization schedule.
  test <- anova(mod,permutations=t(apply(P12,1,order)),model="direct")
  rd[[k]] <- data.frame(predictor=k,R2=mod$CCA$tot.chi/mod$tot.chi,p=test[1,"Pr(>F)"])
}
out(do.call(rbind,reg),"Table_2_spatial_regressions");out(do.call(rbind,rd),"site_RDA")
source("figures.R", encoding="UTF-8")
capture.output(sessionInfo(),file="results/sessionInfo.txt")
stopifnot(abs(sens$primary["site_R2"]-.63475658)<1e-6,abs(sens$primary["PC1"]-.38443626)<1e-6)
# Keep tracked Word markup and replace its three active figures with this run.
if (!requireNamespace("xml2",quietly=TRUE) || !requireNamespace("zip",quietly=TRUE))
 stop("Install xml2 and zip to export the Word manuscript.")
tmp <- tempfile(); dir.create(tmp)
unzip("manuscript/Curitiba_R_review.docx",exdir=tmp)
for(i in 1:3) stopifnot(file.copy(sprintf("results/Figure_%d.png",i),
 file.path(tmp,"word","media",sprintf("final%d.png",i)),overwrite=TRUE))
part <- file.path(tmp,"word","document.xml")
doc <- xml2::read_xml(part)
nodes <- xml2::xml_find_all(doc,"//w:t",xml2::xml_ns(doc))
swap <- function(old,new) {
  for(node in nodes[grepl(old,xml2::xml_text(nodes),fixed=TRUE)])
   xml2::xml_set_text(node,sub(old,new,xml2::xml_text(node),fixed=TRUE))
}
for(pair in list(c("urban_percent","0.658"),c("forest_percent","0.494"),c("distance_km","0.361")))
 swap(pair[2],sprintf("%.3f",rd[[pair[1]]]$p))
dispersion <- anova(betadisper(dist(z),site,type="median"))
swap("0.426)",sprintf("%.3f)",dispersion[1,"Pr(>F)"]))
swap("P = 0.0001).",sprintf("P = %.4f).",sens[["omit Cr+Zn"]]["p"]))
# Update the signed slopes in Table 2 to use the same axis orientation as Figures 2 and 3.
table2 <- tail(xml2::xml_find_all(doc,"//w:tbl",xml2::xml_ns(doc)),1)
rows <- xml2::xml_find_all(table2,"./w:tr",xml2::xml_ns(doc))
tab <- do.call(rbind,reg)
for(i in seq_len(nrow(tab))) {
  cell <- xml2::xml_find_first(rows[i+1],"./w:tc[3]//w:t",xml2::xml_ns(doc))
  xml2::xml_set_text(cell,sprintf("%.4f",tab$slope[i]))
}
xml2::write_xml(doc,part)
zip::zip("results/Curitiba_R_manuscript_tracked.docx",files=list.files(tmp,all.files=TRUE,no..=TRUE),root=tmp,mode="mirror")
