

setwd("E:\\Workspace\\Sand burner\\SandBox_4000s_Data")
setwd("E:\\Workspace\\Sand burner\\4000s_21Feb2019")
setwd("E:\\Workspace\\Sand burner\\032519 data")
setwd("G:\\Workspace\\bigsandy0919")

flist <- Sys.glob("*.lvm")

nfiles <- length(flist)
library(foreach)
library(doParallel)




library(EMD)

library(parallel)


i <- 445



get.sub.mat <- function(i)
{
  library(EMD)
  setwd("G:\\Workspace\\bigsandy0919")
  
  
  
  curfile <- flist[i]
  
  
  ###Get input data 
  
  fsplit <- strsplit(curfile, "_")
  bnum <- fsplit[[1]][7]
  bnum <-  gsub("[.lvm]", "", bnum )
  
  ###Get slope
  
  slope.str <- fsplit[[1]][3]
  sl.rep <- gsub("[deg]", "", slope.str)
  sl.val <- as.numeric(sl.rep)
  
  
  ###Get slpm
  
  slope.str <- fsplit[[1]][5]
  sl.rep <- gsub("[slpm]", "", slope.str)
  # slpm.split  <- strsplit(sl.rep, "x")
  slpm.val <- as.numeric(sl.rep )
  #slpm.val <- as.numeric(sl.rep )
  
  
  ###Get nburners  
  w.str <- fsplit[[1]][4]
  w.str <- gsub("[B]", "",w.str )
  
  w.val <-   w.str 
  
  sb.dat <- read.table(curfile, skip = 24,header = F)
  
  #sb.class <-slpm.split[[1]][1]
  fzd.sub <-  w.str
  
  #sb.class.unique <- unique(sb.class)
  
  tseq1 <- sb.dat[,1]
  
  inc <- 1
  
  
  ###Get burner sequence
  w.str <- fsplit[[1]][6]
  
  burnerseq <- w.str
  
  
  
  ###Section 1
  #nc <- dim(sb.dat)[2]-2
  nc <- 65
  tseq <- sb.dat[,1]
  sb.sub <- sb.dat
  colseq <- 2:nc
  timevec <- tseq
  #
  # #inc <- 1
  # for(curcol in colseq)
  # {
  #    fout <- paste("Plot-", inc, ".png")
  #    png(fout)
  #    plot(	sb.sub[,curcol], type = 'l')
  #    dev.off()
  #    inc <- inc + 1
  #  }
  
  lc.freq <- rep(NA, length(colseq)+1)
  maxspec.freq <- rep(NA, length(colseq))
  cordist.mat <- matrix(0, nrow = length(colseq), ncol = length(colseq))
  inc <- 1
  avg.vec1 <- rep(NA, length(colseq)+1)
  sd.vec1 <- rep(NA, length(colseq)+1)
  
  p550.vec1 <- rep(NA, length(colseq)+1)
  cor.vec <-rep(NA, length(colseq)+1)
  
  for(curcol in colseq)
  {	  cur.sb <- sb.sub[,curcol]
  cur.avg <- mean(cur.sb)
  #using standard deviation
  #cur.sd <- sd(cur.sb)
  #using MAD
  cur.sd <- mad(cur.sb)
  if(cur.sd>0)
  {	  avg.vec1[curcol] <- cur.avg
  sd.vec1[curcol] <- cur.sd
  p550.vec1[curcol] <- sum(cur.sb > 550) / length(cur.sb)
  cur.extrema <- extrema(cur.sb-cur.avg)
  cur.ncross <- cur.extrema$ncross
  intlength <- max(timevec) - min(timevec)
  cur.lc.freq <- cur.ncross / intlength
  lc.freq[inc] <- cur.lc.freq / 2
  ###Level crossing
  ###Spectral density 
  cur.cent <- cur.sb - cur.avg
  cent.ts <- ts(cur.cent, freq = 500)
  cur.spec <- spectrum(cent.ts, plot = F, spans = c(15,5))
  #plot((cur.spec$freq), (cur.spec$spec),type = "l", xlim = c(0, 25))
  max.freq <- cur.spec$freq[which.max(cur.spec$spec)]
  maxspec.freq [inc] <- max.freq
  ###Correlation distance
  inc2 <- 1
  for(tail.col in colseq)
  {	  headvec <- cur.sb 
  tailvec <- sb.sub[,tail.col]
  cordist.mat[inc, inc2] <- cor(headvec, tailvec)
  inc2 <- inc2 + 1
  }
  inc <- inc + 1
  }
  }
  corvec <- cordist.mat[,3]
  ncor <- 1:length(corvec)
  dvec <- (colseq-2)[-1]
  avec <- avg.vec1[ncor]
  svec <- sd.vec1[ncor]
  lvec <- lc.freq[ncor]
  navec <- is.na(avec)
  valvec <- as.logical(1-navec)
  cvec <- corvec[-64]
  avec <- avec[valvec]
  lvec <- rep(0, length(cvec))
  svec <- svec[valvec]
  dseq1 <- seq(4, 152, by = 4) 
  dseq2 <- seq(160, 352, by = 8)
  dseq <- c(dseq1, dseq2)
  
  #temp.list[[i]] <- avg.vec1
  #freq.list[[i]] <- lc.freq
  #corr.list[[i]] <- corvec
  slvec  <- rep(sl.val, length(cvec ))
  evec <- rep(slpm.val, length(cvec ))
  #wsvec<- rep(wind.val, length(cvec ))
  fzvec <- rep(fzd.sub, length(cvec ))
  #widthvec <- rep(w.val, length(cvec ))
  bseqvec <- rep(burnerseq, length(cvec))
  dvec1 <- 2:64
  dvec2 <- na.exclude(dvec1[valvec])
  p550.val <- p550.vec1[dvec1 ]
  cvec <- corvec[dvec]
  outmat <- cbind(avec,lvec , cvec, dseq, slvec, evec,fzvec, p550.val, svec, bseqvec)
  bnumvec <- rep(bnum, dim(outmat)[1])
  rseq <- 1:dim(outmat)[1]
  rvec <- paste(bnumvec, rseq, sep = "-")
  rownames(outmat)<-rvec
  #plot(colseq, corvec)
  #abline(cor.regression[2],cor.regression[1] )
  #lines(colseq, cor.pred)
  #setwd("E:\\Workspace\\Sand burner\\Output-4000s")
  #write.table(lc.freq, "lc-freq-1.txt", row.names = F)
  #write.table(cordist.mat, "cor-dist-1.txt", row.names = F)
  fout <- paste(bnum, "-mad-out.txt", sep="")
  
  
  ###Removing bad TCs
  outmat <- outmat[-c(7, 27, 29),]
  write.table(outmat, fout)  
  
  
  fout <- paste(bnum, "-p550.txt", sep="")
  #write.table(p550.vec1, fout)
  
  
  ###Get residuals, and ACFs
  #tempx <- log(dvec)
  #tempy <- log(avec)
  #templm <- lm(tempy ~ tempx)
  #abline(coef(templm)[1], coef(templm)[2])
  #temp.sp <- smooth.spline()
  print(i / nfiles)
  
  
  dseq <- seq(4, 324, by = 10)
  pseq <- 65:97
  ptseq <- 97:129	 
  
  pmat <- sb.dat[,pseq]
  ptmat <- sb.dat[,ptseq]
  
  pvec <- apply(pmat, 2, "mean")
  ptvec <- apply(ptmat, 2, "mean")
  
  poutmat <- cbind(dseq, pvec, ptvec)
  indseq <- 1:dim(poutmat)[1]
  isvalid <- as.logical(pvec + ptvec < 1000)
  indval <- indseq[isvalid]
  poutmat <- poutmat[indval,]
  fout <- paste(bnum, "-pressure-temp.txt", sep = "")
  write.table(poutmat, fout)
  
  
  
}


cl <- makeCluster(64)
registerDoParallel(cl)


a <- foreach(m = 1:147, .packages = "EMD") %dopar%   get.sub.mat(m)
closeAllConnections()


matflist <- Sys.glob("*.txt")

pflist <- grep("pressure", matflist)
pflist <- matflist[pflist]

plist <- vector("list", length(pflist))

for(i in 1:length(pflist))
{
  curf <- pflist[i]
  curstr <- strsplit(curf, "-")
  curlab <- curstr[[1]][1]
  curmat <- read.table(pflist[i])
  rvec <- 1:dim(curmat)[1]
  rlab <- paste(rvec , curlab, sep = "-")
  rownames(curmat) <- rlab
  plist[[i]] <- curmat
}

pmat <- do.call("rbind", plist)

write.table(pmat, "press-temp-mat.csv", sep = ",")



mflist <- grep("mad", matflist)
mflist <- matflist[mflist]

mlist <- vector("list", length(mflist))

for(i in 1:length(mflist))
{
  curf <- mflist[i]
  curstr <- strsplit(curf, "-")
  curlab <- curstr[[1]][1]
  curmat <- read.table(mflist[i])
  rvec <- 1:dim(curmat)[1]
  rlab <- paste(rvec , curlab, sep = "-")
  rownames(curmat) <- rlab
  mlist[[i]] <- curmat
}

mmat <- do.call("rbind", mlist)

write.table(mmat, "temp-mat.csv", sep = ",")





for(i in 1:70)
{
  get.sub.mat(i)
  
}

flist.in <- flist
firenum.vec <- rep("", nfiles)
for(i in 1:nfiles)
{ #setwd("E:\\Workspace\\Sand burner\\SandBox_4000s_Data")
  setwd("E:\\Workspace\\Sand burner\\4000s_21Feb2019")
  
  curfile <- flist[i]
  
  
  ###Get input data 
  
  fsplit <- strsplit(curfile, "_")
  ftemp <- fsplit[[1]][6]
  ftemp <- gsub("[AA]", "", ftemp)
  ftemp <- gsub("[.lvm]", "", ftemp)
  firenum.vec[i] <- ftemp
  
}

for(i in 1:nfiles)
{	
  
}


setwd("E:\\Workspace\\Sand burner\\Output-4000s")
setwd("E:\\Workspace\\Sand burner\\Output-4000s\\Output-4000-031419")

files.out <- Sys.glob("*.txt")
nfiles.out <- length(files.out)
firenum.vec.out <- rep("", nfiles.out)
for(i in 1:nfiles.out)
{
  curfile <- files.out[i]
  ftemp <- gsub("[-mad-out.txt]", "", curfile)
  ftemp <- gsub("[AA]", "", ftemp )
  firenum.vec.out[i] <- ftemp
}

dvec <-setdiff(as.numeric(firenum.vec),as.numeric(firenum.vec.out))
