
## Person-level results

round(cbind(output_P1$result[,c(2,3,4,7)], output_P1$ratio.gw), 3)

round(cbind(output_P3$result[,c(2,3,4,7)], output_P3$ratio.gw), 3)

round(cbind(output_P4$result[,c(2,3,4,7)], output_P4$ratio.gw), 3)

round(cbind(output_P5$result[,c(2,3,4,7)], output_P5$ratio.gw), 3)

round(cbind(output_P6$result[,c(2,3,4,7)], output_P6$ratio.gw), 3)

round(cbind(output_P7$result[,c(2,3,4,7)], output_P7$ratio.gw), 3)

round(cbind(output_P9$result[,c(2,3,4,7)], output_P9$ratio.gw), 3)

round(cbind(output_P10$result[,c(2,3,4,7)], output_P10$ratio.gw), 3)

round(cbind(output_P11$result[,c(2,3,4,7)], output_P11$ratio.gw), 3)

output_P1$result[which(output_P1$result[,4] == min(output_P1$result[,4])),1]
output_P3$result[which(output_P3$result[,4] == min(output_P3$result[,4])),1]
output_P4$result[which(output_P4$result[,4] == min(output_P4$result[,4])),1]
output_P5$result[which(output_P5$result[,4] == min(output_P5$result[,4])),1]
output_P6$result[which(output_P6$result[,4] == min(output_P6$result[,4])),1]
output_P7$result[which(output_P7$result[,4] == min(output_P7$result[,4])),1]
output_P9$result[which(output_P9$result[,4] == min(output_P9$result[,4])),1]
output_P10$result[which(output_P10$result[,4] == min(output_P10$result[,4])),1]
output_P11$result[which(output_P11$result[,4] == min(output_P11$result[,4])),1]


#Save data
P_names <- paste0("P", 1:11)
#P_names <- matrix(P_names[-c(2,8)], ncol = 1)
P_names <- matrix(P_names, ncol = 1)
colnames(P_names) <- "Participant"
#P_names
mydata <- rbind(
cbind.data.frame(cbind(P_names[1], output_P1$result[,1]), round(cbind(output_P1$result[,c(2,3,4,7)], output_P1$ratio.gw), 3)),
cbind.data.frame(cbind(P_names[3], output_P3$result[,1]), round(cbind(output_P3$result[,c(2,3,4,7)], output_P3$ratio.gw), 3)),
cbind.data.frame(cbind(P_names[4], output_P4$result[,1]), round(cbind(output_P4$result[,c(2,3,4,7)], output_P4$ratio.gw), 3)),
cbind.data.frame(cbind(P_names[5], output_P5$result[,1]), round(cbind(output_P5$result[,c(2,3,4,7)], output_P5$ratio.gw), 3)),
cbind.data.frame(cbind(P_names[6], output_P6$result[,1]), round(cbind(output_P6$result[,c(2,3,4,7)], output_P6$ratio.gw), 3)),
cbind.data.frame(cbind(P_names[7], output_P7$result[,1]), round(cbind(output_P7$result[,c(2,3,4,7)], output_P7$ratio.gw), 3)),
cbind.data.frame(cbind(P_names[9], output_P9$result[,1]), round(cbind(output_P9$result[,c(2,3,4,7)], output_P9$ratio.gw), 3)),
cbind.data.frame(cbind(P_names[10], output_P10$result[,1]), round(cbind(output_P10$result[,c(2,3,4,7)], output_P10$ratio.gw), 3)),
cbind.data.frame(cbind(P_names[11], output_P11$result[,1]), round(cbind(output_P11$result[,c(2,3,4,7)], output_P11$ratio.gw), 3))
)
mydata
write.table(mydata, "Results.xls", row.names=FALSE, sep="\t", dec=',')
#
mydata_1000 <- mydata
mydata_numb_mx <- mydata[,-(c(1,2))]
mydata_numb_vec <- formatC(t(mydata_numb_mx), format = "f", digits = 3) # t(mydata_numb_mx)
mydata_numb_vec[which(t(mydata_numb_mx) > 1000)] <- ">1000"
mydata_numb_vec
mydata_numb_mx <- matrix(mydata_numb_vec, byrow = T, ncol = dim(mydata_numb_mx)[2])
mydata_numb_mx
mydata_1000[,-(c(1,2))] <- mydata_numb_mx
mydata_1000
write.table(mydata_1000, "Results_1000.xls", row.names=FALSE, sep="\t", dec=',')


######################
