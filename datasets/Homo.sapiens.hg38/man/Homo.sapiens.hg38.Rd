\name{Homo.sapiens.hg38}
\docType{data}
\alias{Homo.sapiens.hg38}
\title{Annotation package that collates several annotation resources.}
\description{This data object was automatically created by Author Name. It
  represents a collection of annotation packages that can be used as a
  single object named Homo.sapiens.hg38.  This object can be used with the
  standard four accessor method for all AnnotationDbi objects.  Namely:
  \code{cols}, \code{keytype}, \code{keys} and \code{select}.  Users are
  encouraged to read the vignette from the \code{OrganismDbi} package for
  more details.}
\usage{library(Homo.sapiens.hg38)}
\keyword{datasets}
\examples{
  Homo.sapiens.hg38
  cls <- columns(Homo.sapiens.hg38)
  cls
  cls <- cls[c(1,19,45)]
  kts <- keytypes(Homo.sapiens.hg38)
  kt <- kts[2]
  kts
  ks <- head(keys(Homo.sapiens.hg38, keytype=kts[2]))
  ks
  res <- select(Homo.sapiens.hg38, keys=ks, columns=cls, keytype=kt)
  head(res)
}
