# Usage: Rscript plot_public_figures.R INPUT_AGGREGATE_DIRECTORY OUTPUT_DIRECTORY
# Aggregate-only scientific figures. R grid and grDevices are the exclusive graphics backend.
# Structural adaptation of the supplied base R results workflow. No patient inputs or simulation.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
input <- normalizePath(args[1], mustWork=TRUE)
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
out <- normalizePath(out, mustWork=TRUE)
library(grid)
stopifnot(capabilities("cairo"),capabilities("png"))
f2 <- read.csv(file.path(input,"figure2-source-data.csv"),stringsAsFactors=FALSE)
fs <- read.csv(file.path(input,"figureS1-source-data.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(f2)==6L,nrow(fs)==16L)
d <- rbind(f2,fs)
stopifnot(sum(duplicated(d))==2L,identical(duplicated(d),duplicated(d[c("center","analysis","version")])))
d <- d[!duplicated(d[c("center","analysis","version")]),,drop=FALSE]
flow <- read.csv(file.path(input,"figure1-source-data.csv"),stringsAsFactors=FALSE)
cs <- read.csv(file.path(input,"figure3-source-data.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(d)==20L,all(d$n00+d$n01+d$n10+d$n11==d$n),
          all(is.finite(d$discordance)),all(d$discordance_lower<=d$discordance_upper),
          nrow(flow)==11L,nrow(cs)==32L)
centers <- c("mover","vitaldb")
center_names <- c(mover="MOVER",vitaldb="VitalDB")
method_order <- c("tail_plus_bins","tail_only","lower_bounded_round_0p2")
method_labels <- c("Quantile coarsening","Tail clipping","0.2 mg/dL grid")
cols <- c(primary="#276B8E",secondary="#78858E",up="#A95C34",down="#276B8E",ink="#20262B",grid="#DDE2E5")
getrow <- function(center,analysis,version) {
  z <- d[d$center==center & d$analysis==analysis & d$version==version,,drop=FALSE]
  stopifnot(nrow(z)==1L); z
}
width_mm <- 170
y_scale <- 1
font_defaults <- gpar(fontfamily="Helvetica",fontsize=8)
ptmin <- 7.2
text_extent <- list()
txt <- function(label,x,y,size=8,face="plain",col=cols["ink"],just="centre",rot=0) {
  stopifnot(size>=ptmin)
  g <- textGrob(label,x=unit(x,"mm"),y=unit(y*y_scale,"mm"),just=just,rot=rot,
                gp=gpar(fontfamily=font_defaults$fontfamily,fontsize=size,fontface=face,col=col,lineheight=1.05))
  grid.draw(g)
}
ln <- function(x0,y0,x1,y1,col=cols["ink"],lwd=.65,lty=1,arrow=FALSE) {
  grid.segments(unit(x0,"mm"),unit(y0*y_scale,"mm"),unit(x1,"mm"),unit(y1*y_scale,"mm"),
                gp=gpar(col=col,lwd=lwd,lty=lty),
                arrow=if(arrow)grid::arrow(length=unit(1.1,"mm"),type="closed") else NULL)
}
box <- function(x,y,w,h,fill="white",border="#9BA5AB",lwd=.65,lty=1) {
  grid.rect(unit(x,"mm"),unit(y*y_scale,"mm"),unit(w,"mm"),unit(h*y_scale,"mm"),
            gp=gpar(fill=fill,col=border,lwd=lwd,lty=lty))
}
mark <- function(x,y,pch=19,col=cols["primary"],fill=col,size=1.8) {
  grid.points(unit(x,"mm"),unit(y*y_scale,"mm"),pch=pch,size=unit(size,"mm"),gp=gpar(col=col,fill=fill,lwd=.9))
}
axis_x <- function(x0,x1,y,limits,at,size=7.5) {
  X <- function(v)x0+(v-limits[1])/diff(limits)*(x1-x0)
  ln(x0,y,x1,y,lwd=.8)
  for(v in at){ln(X(v),y,X(v),y-1.1,lwd=.65);txt(format(v,trim=TRUE),X(v),y-3.4,size=size)}
  invisible(X)
}
interval <- function(X,lo,mid,hi,y,col=cols["primary"],pch=19) {
  stopifnot(is.finite(lo),is.finite(mid),is.finite(hi),lo<=hi)
  ln(X(lo),y,X(hi),y,col=col,lwd=1.05)
  ln(X(lo),y-.7,X(lo),y+.7,col=col,lwd=.7)
  ln(X(hi),y-.7,X(hi),y+.7,col=col,lwd=.7)
  mark(X(mid),y,pch=pch,col=col,fill=if(pch==21)"white" else col)
}
# Every PDF text glyph is at least 7.2 pt, and source values stay unrounded until labeling.
draw_flow <- function() {
  grid.newpage()
  txt("a  Cohort selection",4,193,size=10,face="bold",just="left")
  txt("Retained cases / patients",4,186,size=7.5,just="left")
  txt("MOVER",101,186,size=7.5,face="bold");txt("VitalDB",124,186,size=7.5,face="bold")
  txt("Excluded at this step",150,190,size=7.3)
  txt("MOVER",142,185.5,size=7.3);txt("VitalDB",162,185.5,size=7.3)
  labels <- c("Unique released operation cases","Nonconflicting identity cases","One operation in complete source",
              "Adult and general anesthesia","Valid anesthesia times","Valid admission / discharge ordering",
              "Eligible noncardiac procedure mapping","No renal-status screen exclusion","Preoperative and postoperative SCr coverage")
  yy <- seq(177,97,by=-10)
  for(i in seq_len(9L)) {
    box(66.5,yy[i],125,7.2,fill=if(i==9)"#F0F4F6" else "white")
    txt(labels[i],6.2,yy[i],size=7.3,just="left")
    txt(as.character(flow$mover_remaining[i]),101,yy[i],size=8)
    txt(as.character(flow$vitaldb_remaining[i]),124,yy[i],size=8)
    if(i>1L){
      ln(66.5,yy[i-1]-3.6,66.5,yy[i]+3.8,arrow=TRUE)
      txt(as.character(flow$mover_excluded[i]),142,yy[i],size=7.6)
      txt(as.character(flow$vitaldb_excluded[i]),162,yy[i],size=7.6)
    }
  }
  ln(66.5,93.4,66.5,90);ln(36,90,126,90);ln(36,90,36,85.5,arrow=TRUE);ln(126,90,126,85.5,arrow=TRUE)
  box(39,78,70,14,fill="#DFEAF0",border=cols["primary"],lwd=1)
  txt("Primary E48 cohorts",39,82,size=8,face="bold")
  txt(sprintf("MOVER %d  |  VitalDB %d",flow$mover_remaining[10],flow$vitaldb_remaining[10]),39,77,size=8)
  txt("At least one target with a prior within 48 h",39,72.8,size=7.2)
  box(127,78,77,14,fill="white",border="#7B8790",lty=2)
  txt("Parallel E7 ratio-only cohorts",127,82,size=8,face="bold")
  txt(sprintf("MOVER %d  |  VitalDB %d",flow$mover_remaining[11],flow$vitaldb_remaining[11]),127,77,size=8)
  txt("Separate branch; prior within 7 days",127,72.8,size=7.2)
  txt(sprintf("No E48 target: %d / %d",flow$mover_excluded[10],flow$vitaldb_excluded[10]),39,66.5,size=7.5)
  txt(sprintf("No E7 target: %d / %d",flow$mover_excluded[11],flow$vitaldb_excluded[11]),127,66.5,size=7.5)
  txt("b  Paired numerical comparison",4,58,size=10,face="bold",just="left")
  box(85,47,162,10,fill="#F5F6F7")
  txt("Primary patients, retained observation times and prior-window indices held fixed",85,49,size=8)
  txt("Each center fits its own mappings to its retained reference measurements",85,45,size=7.5)
  ln(85,42,85,39);ln(42,39,127,39);ln(42,39,42,35.5,arrow=TRUE);ln(127,39,127,35.5,arrow=TRUE)
  box(42,29,76,12,fill="white")
  txt("Original released values (I)",42,30.5,size=8.3,face="bold")
  txt("Computational reference",42,26,size=7.5)
  box(127,29,77,12,fill="#DFEAF0",border=cols["primary"])
  txt("Quantile coarsening (B)",127,30.5,size=8.3,face="bold")
  txt("Tail clipping (T) and grid rounding (G) secondary",127,26,size=7.2)
  ln(42,23,42,20);ln(127,23,127,20);ln(42,20,127,20);ln(85,20,85,17,arrow=TRUE)
  box(85,10.5,162,11,fill="#F5F6F7")
  txt("Primary I versus B paired labels: n00, n01, n10, n11",85,12.5,size=8.1,face="bold")
  txt("Discordance D = (n01 + n10) / N; report both directions and net change separately",85,8,size=7.8)
  txt("Flow counts start at unique cases. Later eligible cohorts contain one operation per patient.",4,1.8,size=7.2,just="left")
}
draw_main <- function() {
  grid.newpage()
  txt("a  Paired label discordance",4,113,size=9.5,face="bold",just="left")
  txt("b  Direction of primary changes",94,113,size=9.5,face="bold",just="left")
  X <- axis_x(41,90,16,c(0,14.8),seq(0,14,2))
  for(v in seq(0,14,2))ln(X(v),17,X(v),99,col=cols["grid"],lwd=.5)
  for(j in seq_along(centers)){
    z0 <- getrow(centers[j],"main","tail_plus_bins")
    yy <- if(j==1)c(88,77,66) else c(48,37,26)
    txt(sprintf("%s  n = %d",center_names[centers[j]],z0$n),41,max(yy)+10,size=8.4,face="bold",just="left")
    for(k in seq_along(method_order)){
      z <- getrow(centers[j],"main",method_order[k]);cl <- if(k==1)cols["primary"] else cols["secondary"]
      txt(method_labels[k],3.5,yy[k],size=7.7,just="left")
      interval(X,z$discordance_lower*100,z$discordance*100,z$discordance_upper*100,yy[k],cl,if(k==1)19 else 21)
    }
  }
  txt("Patients with changed label (%)",65.5,7,size=7.6)
  XP <- axis_x(115,166,16,c(-4.1,4.1),seq(-4,4,2))
  ln(XP(0),17,XP(0),99,col="#57636C",lwd=.8,lty=2)
  for(j in seq_along(centers)){
    z <- getrow(centers[j],"main","tail_plus_bins"); yy <- if(j==1)c(87,70) else c(47,30)
    txt(center_names[centers[j]],105,max(yy)+11,size=8.5,face="bold",just="left")
    txt("1 to 0",112,yy[1],size=7.7,just="right");txt("0 to 1",112,yy[2],size=7.7,just="right")
    interval(XP,-z$downward_upper*100,-z$downward*100,-z$downward_lower*100,yy[1],cols["down"],17)
    interval(XP,z$upward_lower*100,z$upward*100,z$upward_upper*100,yy[2],cols["up"],19)
    txt(sprintf("%d/%d (%.2f%%)",z$n10,z$n,z$downward*100),XP(-z$downward*100),yy[1]-5.7,size=7.3)
    txt(sprintf("%d/%d (%.2f%%)",z$n01,z$n,z$upward*100),XP(z$upward*100),yy[2]-5.7,size=7.3)
  }
  txt("Patients with directional change (%)",140.5,7,size=7.6)
  txt("Lines show 2.5th-97.5th percentile ranges from 5000 patient resamples; population coverage is unestablished.",85,1.8,size=7.2)
}
states <- c("neither","ratio_only","absolute_only","both")
state_labels <- c("Neither","Ratio\nonly","Absolute\nonly","Both")
cs$original_index <- match(cs$original_state,states)
cs$fixed_B_index <- match(cs$fixed_B_state,states)
cs$composite_flip <- (cs$original_index==1L)!=(cs$fixed_B_index==1L)
cs$state_change <- cs$original_index!=cs$fixed_B_index
cs$transition_class <- ifelse(!cs$state_change,"Unchanged state",ifelse(cs$composite_flip,"Composite flip","Positive-state change"))
for(center in centers){
 z <- cs[cs$dataset==center,]; ref <- getrow(center,"main","tail_plus_bins")
 stopifnot(sum(z$patients)==ref$n, sum(z$patients[z$composite_flip])==ref$discordance_count)
 stopifnot(sum(z$patients[z$state_change])==if(center=="mover")225 else 106)
}
draw_states <- function() {
 grid.newpage()
 for(j in seq_along(centers)) {
  center <- centers[j]; z <- cs[cs$dataset==center,];ref <- getrow(center,"main","tail_plus_bins")
  left <- if(j==1)27 else 112;cell <- 12.4; top <- 85
  txt(sprintf("%s  %s  n = %d",letters[j],center_names[center],ref$n),left-23,107,size=9.5,face="bold",just="left")
  txt("After fixed quantile coarsening",left+2*cell,99,size=7.8)
  for(k in seq_len(4)){
   txt(state_labels[k],left+(k-.5)*cell,91.5,size=7.4)
   txt(state_labels[k],left-2.2,top-(k-.5)*cell,size=7.4,just="right")
  }
  for(i in seq_len(nrow(z))){
   zz <- z[i,]; cx <- left+(zz$fixed_B_index-.5)*cell;cy <- top-(zz$original_index-.5)*cell
   flip <- zz$composite_flip && zz$patients>0; pos <- zz$state_change && !zz$composite_flip && zz$patients>0
   fill <- if(flip)"#DFEAF0" else if(pos)"#F1E8DF" else if(!zz$state_change)"#F1F3F4" else "white"
   border <- if(flip)cols["primary"] else if(pos)"#80664C" else "#C9D0D5"
   box(cx,cy,cell-.8,cell-.8,fill=fill,border=border,lwd=if(flip)1.05 else .6,lty=if(pos)2 else 1)
   txt(as.character(zz$patients),cx,cy,size=9,face=if(flip)"bold" else "plain",col=if(zz$patients==0)"#7B848A" else cols["ink"])
  }
  flipn <- sum(z$patients[z$composite_flip]);staten <- sum(z$patients[z$state_change]);positiven <- staten-flipn
  txt(sprintf("Composite flips  %d",flipn),left+2*cell,29.5,size=8.5,face="bold")
  txt(sprintf("Any state change  %d",staten),left+2*cell,24.4,size=8.3)
  txt(sprintf("%d flips + %d positive-state changes",flipn,positiven),left+2*cell,19.3,size=7.5)
 }
 txt("Original states are rows; transformed states are columns. Cell labels are patient counts.",85,12.5,size=7.5)
 box(7,6.2,3.2,3.2,fill="#DFEAF0",border=cols["primary"],lwd=1);txt("Composite flip",10,6.2,size=7.2,just="left")
 box(54,6.2,3.2,3.2,fill="#F1E8DF",border="#80664C",lty=2);txt("Positive-state change",57,6.2,size=7.2,just="left")
 box(115,6.2,3.2,3.2,fill="#F1F3F4",border="#C9D0D5");txt("Unchanged state",118,6.2,size=7.2,just="left")
 txt("Post-result exploratory analysis. Both criteria may be met at different E48 times. No new intervals or tests.",85,2.0,size=7.2)
}
modes <- c("main","parallel7","S1_no_diagnosis","S2_target_after5min","S3_exclude_extremes","S4_exclude_over24h","S5_exclude_conflicts","S6_all_postop_E48")
mlabels <- c("Primary","Parallel 7-day ratio only*","S1  No diagnosis screen","S2  Targets after 5 min","S3  Exclude extreme values","S4  Exclude anesthesia >24 h","S5  Exclude time conflicts","S6  All observed points in E48")
draw_sensitivity <- function() {
 grid.newpage()
 txt("Analysis",3,111,size=8.5,face="bold",just="left")
 txt("a  MOVER",74,113,size=9.5,face="bold",just="left")
 txt("b  VitalDB",129,113,size=9.5,face="bold",just="left")
 txt("n",66,104,size=7.7);txt("n",121,104,size=7.7)
 yy <- seq(97,27,by=-10)
 for(k in seq_along(modes))txt(mlabels[k],3,yy[k],size=7.7,just="left")
 for(j in seq_along(centers)){
   x0 <- if(j==1)74 else 129;x1 <- x0+37
   XX <- axis_x(x0,x1,19,c(0,5.6),0:5,size=7.4)
   for(v in 0:5)ln(XX(v),20,XX(v),102,col=cols["grid"],lwd=.5)
   main <- getrow(centers[j],"main","tail_plus_bins")
   ln(XX(main$discordance*100),20,XX(main$discordance*100),102,col="#A1B7C3",lty=3,lwd=.8)
   for(k in seq_along(modes)){
     z <- getrow(centers[j],modes[k],"tail_plus_bins")
     txt(as.character(z$n),x0-4,yy[k],size=7.6,just="right")
     interval(XX,z$discordance_lower*100,z$discordance*100,z$discordance_upper*100,yy[k],
              if(k==1)cols["primary"] else cols["secondary"],if(k==1)19 else 21)
   }
 }
 txt("Patients with changed label (%)",121,10.5,size=7.7)
 txt("Quantile coarsening. Lines show resampling ranges. Each analysis refits its mapping.",85,6.2,size=7.2)
 txt("*Parallel E7 changes population, rule and reference pool. S6 does not establish complete seven-day surveillance.",85,2.1,size=7.2)
}
# Native Cairo PDF preserves selectable text and embeds fonts. Raster outputs use physical mm dimensions.
export <- function(stem,fun,w_mm,h_mm) {
  w <- w_mm/25.4;h <- h_mm/25.4
  y_scale <<- if(stem=="figure1-cohort-and-paired-design")185/200 else 1
  if(requireNamespace("svglite",quietly=TRUE)) {
    svglite::svglite(file.path(out,paste0(stem,".svg")),width=w,height=h,
                    system_fonts=list(sans="Helvetica"))
    fun();dev.off()
  }
  grDevices::cairo_pdf(file.path(out,paste0(stem,".pdf")),width=w,height=h,family="Helvetica",pointsize=8)
  fun();dev.off()
  grDevices::png(file.path(out,paste0(stem,".png")),width=w_mm,height=h_mm,units="mm",res=300,type="cairo",pointsize=8)
  fun();dev.off()
  tiff_path <- file.path(out,paste0(stem,".tiff"))
  if(capabilities("tiff")) {
    grDevices::tiff(tiff_path,width=w_mm,height=h_mm,units="mm",res=600,type="cairo",compression="lzw",pointsize=8)
    fun();dev.off()
  } else if(capabilities("aqua") && nzchar(Sys.which("sips"))) {
    grDevices::tiff(tiff_path,width=w_mm,height=h_mm,units="mm",res=600,type="quartz",pointsize=8)
    fun();dev.off()
    # R draws the TIFF natively. System ImageIO only losslessly encodes LZW and restores density metadata.
    status <- system2(Sys.which("sips"),c("--setProperty","formatOptions","lzw",
                         "--setProperty","dpiWidth","600","--setProperty","dpiHeight","600",shQuote(tiff_path)),stdout=FALSE)
    stopifnot(status==0L)
  } else warning("Optional TIFF omitted: no native TIFF or Quartz/ImageIO export route")
}
export("figure1-cohort-and-paired-design",draw_flow,width_mm,185)
export("figure2-paired-label-changes",draw_main,width_mm,120)
export("figure3-criterion-state-transitions",draw_states,width_mm,114)
export("figureS1-prespecified-analyses",draw_sensitivity,width_mm,120)
writeLines(c(R.version.string,"Backend: R grid/grDevices only","PDF: Cairo embedded selectable text","PNG: 300 dpi; TIFF: 600 dpi LZW","Optional SVG omitted because svglite is unavailable; native Cairo SVG outlines glyphs.","All input observations are frozen aggregates. No patient calculations or resampling."),file.path(out,"render-environment.txt"))
