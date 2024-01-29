library("circlize")
setwd("C:/Users/All users.DESKTOP-B09KNMV/Documents/FYP_lynn/metagenome analysis/chordplot")
EV<-read.csv("EV_csv.csv",row.names = "X")
EV<-as.matrix(EV)

MF<-read.csv("MF_csv.csv",row.names = "X")
MF<-as.matrix(MF)


PB<-read.csv("PB_csv.csv",row.names = "X")
PB<-as.matrix(PB)



colors <- c(binding = "#848ccf", cellulose = "#93b5e1",
            hemicellulose = "#394989", lignin = "#8fcfd1",
            pectin = "#436f8a", Others = "#F7C530",
            Unclassified = "#F0ECE3", Proteobacteria = "#FE9801",
            Actinobacteria = "#D89CF6", Bacteroidetes="#CCDA46", 
            Planctomycetes="#2C786C",Tenericutes="#7D5A5A")

par(mfrow = c(1, 3))
par(cex = 1.3, mar = c(0, 0, 0, 0))

chordDiagram(EV,grid.col = colors)
chordDiagram(MF,grid.col = colors)
chordDiagram(PB,grid.col = colors)

circos.clear()
