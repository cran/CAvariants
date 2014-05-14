cabasicresults<-setClass("cabasicresults",
representation(
  RX="matrix", CX="matrix", Rweights="matrix", Cweights="matrix",
  Raxes="matrix", Caxes="matrix", r="numeric", mu="numeric",mu2="numeric",catype="character",
tau="numeric",tauDen="numeric",Z="matrix",ZtZ="matrix",tZZ="matrix"))
