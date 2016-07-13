libraries = c("tuneR")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

libraries = c("audio")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

library(tuneR)
library(audio)
w <-noise(kind = c("white"))
p <-noise(kind = c("pink"))
b <-noise(kind=c("power"))
par(mfrow=c(3,1))
plot(w,main="white noise")
plot(p,main="pink noise")
plot(p,main="blue noise")
writeWave(w,"w.wav")#writes pink noise on your hard drive
writeWave(p,"p.wav")#loads `audio` package to use `load.wave` function
writeWave(b,"b.wav")