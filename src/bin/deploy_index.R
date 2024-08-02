library(distill)
library(connectapi)

args <- commandArgs(trailingOnly=TRUE)
site_dir <- args[1]
site_name <- args[2]

setwd(site_dir)

rsconnect::addConnectServer(url = "https://connect.sparkds.io/", name = "connect", quiet = TRUE)
rsconnect::connectApiUser(account = "Chao.Di@sparktx.com", server = "connect", apiKey = "za4RPQdH58QuGHPTZbHvU087nF4TlkD3")
rsconnect::deploySite(siteName = site_name, account = "Chao.Di@sparktx.com", server = "connect", launch.browser = FALSE, render = "server")
