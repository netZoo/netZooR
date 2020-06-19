library(covr)
fileName <- '~/token.txt'
token=readChar(fileName, file.info(fileName)$size);token=gsub("[\r\n]", "", token)
codecov(token=token)
