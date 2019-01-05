
rm(list=ls())
library(data.table)
library(ggplot2)
library(scales)



# assume Parallel Trends but violate No Anticipation

dd <- data.table(year = 1990:2005)

dd[year < 1995, cohort1995noantic := .7]
dd[year >= 1995, cohort1995noantic := .6]
dd[, cohort1995noantic_never := .7]

dd[year < 1993, cohort1995antic := .7]
dd[year >= 1993 & year <= 1994, cohort1995antic := .65]
dd[year >= 1995, cohort1995antic := .6]

dd[year < 1995+5, cohort2000noantic := .75]
dd[year >= 1995+5, cohort2000noantic := .65]
dd[, cohort2000noantic_never := .75]

dd[year < 1993+5, cohort2000antic := .75]
dd[year >= 1993+5 & year <= 1994+5, cohort2000antic := .7]
dd[year >= 1995+5, cohort2000antic := .65]

for(var in c("cohort1995noantic","cohort1995antic","cohort2000noantic","cohort2000antic","cohort1995noantic_never","cohort2000noantic_never")){
  dd[, (var) := get(var) - 0.005*(year-1990)]
}


gg <- ggplot(aes(x=year),data=dd) +
  geom_line(aes(y=cohort1995noantic), color = "black") + geom_point(aes(y=cohort1995noantic,colour="E=1995")) +
  geom_line(aes(y=cohort2000noantic), color = "red") + geom_point(aes(y=cohort2000noantic,colour="E=2000")) +
  geom_line(aes(y=cohort1995noantic_never), color = "black", linetype="dashed") +
  geom_line(aes(y=cohort2000noantic_never), color = "red", linetype="dashed") +
  scale_y_continuous(breaks= pretty_breaks(),limits=c(0.4,0.8)) +
  scale_x_continuous(breaks= pretty_breaks(),limits=c(1990,2005)) +
  theme_bw(base_size = 16) +
  scale_colour_manual("Husband Died", breaks = c("E=1995", "E=2000"), values = c("black", "red")) +
  annotate(geom="text",x=1994,y=.56,label="Treatment") +
  annotate(geom="text",x=1999,y=.59,label="Treatment",color="red") +
  annotate(geom="text",x=1994.5,y=.44,label="E=2000 is valid control group\nfor E=1995 in shaded area",colour="blue") +
  annotate("rect", xmin = 1990, xmax = 1999, ymin = .4, ymax = .8, alpha = .2) +
  ylab("Wife's Mean Labor Supply") + xlab("Year")
print(gg)
ggsave(gg,file="fadlon_ideal.pdf",width=8,height=5)


gg <- ggplot(aes(x=year),data=dd) +
  geom_line(aes(y=cohort1995antic), color = "black") + geom_point(aes(y=cohort1995antic,colour="E=1995")) +
  geom_line(aes(y=cohort2000antic), color = "red") + geom_point(aes(y=cohort2000antic,colour="E=2000")) +
  geom_line(aes(y=cohort1995noantic_never), color = "black", linetype="dashed") +
  geom_line(aes(y=cohort2000noantic_never), color = "red", linetype="dashed") +
  scale_y_continuous(breaks= pretty_breaks(),limits=c(0.4,0.8)) +
  scale_x_continuous(breaks= pretty_breaks(),limits=c(1990,2005)) +
  theme_bw(base_size = 16) +
  scale_colour_manual("Husband Died", breaks = c("E=1995", "E=2000"), values = c("black", "red")) +
  annotate(geom="text",x=1992,y=.61,label="Anticipation") +
  annotate(geom="text",x=1997,y=.65,label="Anticipation",color="red") +
  annotate(geom="text",x=1994,y=.55,label="Treatment") +
  annotate(geom="text",x=1999,y=.59,label="Treatment",color="red") +
  annotate(geom="text",x=1993.5,y=.45,label="E=2000 is valid control group\nfor E=1995 in shaded area",colour="blue") +
  annotate("rect", xmin = 1990, xmax = 1997, ymin = .4, ymax = .8, alpha = .2) +
  ylab("Wife's Mean Labor Supply") + xlab("Year")
print(gg)
ggsave(gg,file="fadlon_anticipation.pdf",width=8,height=5)




# assume No Anticipation but violate Parallel Trends

dd <- data.table(year = 1990:2005)
dd[, cohort1995 := .75 + (1990-year)*.02 - 0.1*(year>= 1995)]
dd[, cohort2000 := .7 + (1990-year)*.01 - 0.1*(year>= 2000)]
dd[, cohort1995resid := .75 + (1990-year)*.02 + - 0.1*(year>= 1995)]
dd[, cohort2000resid := .7 + (1990-year)*.02 + - 0.1*(year>= 2000)]

gg <- ggplot(aes(x=year),data=dd) +
  geom_line(aes(y=cohort1995), color = "black") + geom_point(aes(y=cohort1995,colour="E=1995")) +
  geom_line(aes(y=cohort2000), color = "red") + geom_point(aes(y=cohort2000,colour="E=2000")) +
  scale_y_continuous(breaks= pretty_breaks(),limits=c(0.3,0.8)) +
  scale_x_continuous(breaks= pretty_breaks(),limits=c(1990,2005)) +
  theme_bw(base_size = 16) +
  scale_colour_manual("Husband Died", breaks = c("E=1995", "E=2000"), values = c("black", "red")) +
  annotate(geom="text",x=1993.5,y=.56,label="Treatment") +
  annotate(geom="text",x=2001.5,y=.52,label="Treatment",color="red") +
  ylab("Wife's Mean Labor Supply") + xlab("Year")
print(gg)
ggsave(gg,file="fadlon_lineardiff.pdf",width=8,height=5)


gg <- ggplot(aes(x=year),data=dd) +
  geom_line(aes(y=cohort1995resid), color = "black") + geom_point(aes(y=cohort1995resid,colour="E=1995")) +
  geom_line(aes(y=cohort2000resid), color = "red") + geom_point(aes(y=cohort2000resid,colour="E=2000")) +
  scale_y_continuous(breaks= pretty_breaks(),limits=c(0.29,0.8)) +
  scale_x_continuous(breaks= pretty_breaks(),limits=c(1990,2005)) +
  theme_bw(base_size = 16) +
  scale_colour_manual("Husband Died", breaks = c("E=1995", "E=2000"), values = c("black", "red")) +
  annotate(geom="text",x=1994,y=.53,label="Treatment") +
  annotate(geom="text",x=2000.5,y=.53,label="Treatment",color="red") +
  annotate(geom="text",x=1994.5,y=.45,label="E=2000 is valid control group\nfor E=1995 in shaded area",colour="blue") +
  annotate("rect", xmin = 1990, xmax = 1999, ymin = .29, ymax = .8, alpha = .2) +
  ylab("Wife's Mean Labor Supply") + xlab("Year")
print(gg)
ggsave(gg,file="fadlon_lineardiffresid.pdf",width=8,height=5)



# violate both No Anticipation and Parallel Trends

dd <- data.table(year = 1990:2000)

dd[year < 1992, cohort1995 := .7 - 0.03*(year-1990)]
dd[year==1992, cohort1995 := .7 - 0.03*(year-1990) + 0.005]
dd[year==1993, cohort1995 := .7 - 0.03*(year-1990) + 0.015]
dd[year==1994, cohort1995 := .7 - 0.03*(year-1990) + 0.025]
dd[year >= 1995, cohort1995 := .7 - 0.03*(year-1990) - 0.05]
dd[, cohort1995true := .7 - 0.03*(year-1990)]

dd[year < 1997, cohort2000 := .7 - 0.03*(year-1990) + 0.0015*(year-1990)^2]
dd[year == 1997, cohort2000 := .7 - 0.03*(year-1990) + 0.0015*(year-1990)^2 + 0.005]
dd[year == 1998, cohort2000 := .7 - 0.03*(year-1990) + 0.0015*(year-1990)^2 + 0.015]
dd[year == 1999, cohort2000 := .7 - 0.03*(year-1990) + 0.0015*(year-1990)^2 + 0.025]
dd[year >= 2000, cohort2000 := .7 - 0.03*(year-1990) + 0.0015*(year-1990)^2 - 0.05]
dd[, cohort2000true := .7 - 0.03*(year-1990) + 0.0015*(year-1990)^2]

gg <- ggplot(aes(x=year),data=dd) +
  geom_line(aes(y=cohort1995), color = "black") + geom_point(aes(y=cohort1995,colour="E=1995")) +
  geom_line(aes(y=cohort2000), color = "red") + geom_point(aes(y=cohort2000,colour="E=2000")) +
  geom_line(aes(y=cohort1995true), color = "black", linetype="dashed") +
  geom_line(aes(y=cohort2000true), color = "red", linetype="dashed") +
  scale_y_continuous(breaks= pretty_breaks(),limits=c(0.3,0.71)) +
  scale_x_continuous(breaks= pretty_breaks(),limits=c(1990,2000)) +
  theme_bw(base_size = 16) +
  scale_colour_manual("Husband Died", breaks = c("E=1995", "E=2000", "E=1995 true", "E=2000 true"), values = c("black", "red","black","red")) +
  annotate(geom="text",x=1994,y=.5,label="Treatment") +
  annotate(geom="text",x=1999.5,y=.47,label="Treatment",color="red") +
  ylab("Wife's Mean Labor Supply") + xlab("Year")
print(gg)
ggsave(gg,file="fadlon_neitherholds.pdf",width=8,height=5)


