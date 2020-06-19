library(covr)
fileName <- '/home/ubuntu/token.txt'
token=readChar(fileName, file.info(fileName)$size);token=gsub("[\r\n]", "", token)
codecov(token=token)
