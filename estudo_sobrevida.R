### ESTUDO DE SOBREVIV�NCIA EM 10 ANOS DE
### MULHERES COM C�NCER DE MAMA EM JUIZ DE FORA
### DIAGNOSTICADAS ENTRE 2003-2005

### ANALISES EXPLORAT�RIAS

## RESUMO NUM�RICO DO BANCO
summary(bancomg)

## AN�LISE COM FOCO NA VARI�VEL RA�A/COR (BRANCA/N�O BRANCA)

# Rela��o entre cor e faixa et�ria
cor.fxidade <- table(bancomg$id_cat1,bancomg$cor1)
cor.fxidade
round(prop.table(cor.fxidade,1),3)
round(prop.table(cor.fxidade,2),3)
summary(cor.fxidade)
fisher.test(cor.fxidade)

# Rela��o entre cor e renda
cor.renda <- table(bancomg$rendacat,bancomg$cor1)
cor.renda
round(prop.table(cor.renda,1),3)
round(prop.table(cor.renda,2),3)
summary(cor.renda)
fisher.test(cor.renda)

# Rela��o entre cor e natureza do servi�o
cor.servico <- table(bancomg$servico,bancomg$cor1)
cor.servico
round(prop.table(cor.servico,1),3)
round(prop.table(cor.servico,2),3)
summary(cor.servico)
fisher.test(cor.servico)

# Rela��o entre cor e atraso no in�cio do tratamento
cor.atraso <- table(bancomg$atraso2,bancomg$cor1)
cor.atraso
round(prop.table(cor.atraso,1),3)
round(prop.table(cor.atraso,2),3)
summary(cor.atraso)
fisher.test(cor.atraso)

# Rela��o entre cor e estadiamento do tumor
cor.estadio <- table(bancomg$estadio,bancomg$cor1)
cor.estadio
round(prop.table(cor.estadio,1),3)
round(prop.table(cor.estadio,2),3)
summary(cor.estadio)
fisher.test(cor.estadio)

# Rela��o entre cor e �bito - 10 anos
cor.status10 <- table(bancomg$status10,bancomg$cor1)
cor.status10
round(prop.table(cor.status10,1),3)
round(prop.table(cor.status10,2),3)
summary(cor.status10)
fisher.test(cor.status10)


## AN�LISE COM FOCO NA VARI�VEL RENDA (M�DIA DO SETOR CENSIT�RIO, EM SAL�RIOS M�NIMOS)

# Rela��o entre renda e faixa et�ria
renda.idade <- table(bancomg$id_cat1,bancomg$rendacat)
renda.idade
round(prop.table(renda.idade,1),3)
round(prop.table(renda.idade,2),3)
summary(renda.idade)
fisher.test(renda.idade)

# Rela��o entre renda e renda raca/cor
renda.cor <- table(bancomg$cor1,bancomg$rendacat)
renda.cor
round(prop.table(renda.cor,1),3)
round(prop.table(renda.cor,2),3)
summary(renda.cor)
fisher.test(renda.cor)

# Rela��o entre renda e natureza do servi�o
renda.servico <- table(bancomg$servico,bancomg$rendacat)
renda.servico
round(prop.table(renda.servico,1),3)
round(prop.table(renda.servico,2),3)
summary(renda.servico)
fisher.test(renda.servico)

# Rela��o entre renda e atraso no in�cio do tratamento
renda.atraso <- table(bancomg$atraso2,bancomg$rendacat)
renda.atraso
round(prop.table(renda.atraso,1),3)
round(prop.table(renda.atraso,2),3)
summary(renda.atraso)
fisher.test(renda.atraso)

# Rela��o entre renda e estadiamento do tumor
renda.estadio <- table(bancomg$estadio,bancomg$rendacat)
renda.estadio
round(prop.table(renda.estadio,1),3)
round(prop.table(renda.estadio,2),3)
summary(renda.estadio)
fisher.test(renda.estadio)

# Rela��o entre renda e �bito em 10 anos
renda.status10 <- table(bancomg$status10,bancomg$rendacat)
renda.status10
round(prop.table(renda.status10,1),3)
round(prop.table(renda.status10,2),3)
summary(renda.status10)
fisher.test(renda.status10)

### AN�LISE N�O PARAM�TRICA - KAPLAN-MEIER

library(survival)

# Criando objeto tipo sobreviv�ncia
y <- Surv(bancomg$tempo10,bancomg$status10)

# Estimando a curva de sobreviv�ncia pelo Kaplan-Meier
KM <- survfit(y ~ 1, data=bancomg)
print(KM, print.rmean=TRUE)
options(max.print=10000)
summary(KM)
plot(KM, xlab = "Tempo (dias)", ylab = "S(t)", main = "Geral", mark.time = F)

## Teste Log-rank (argumento rho=0 � o default, n�o precisa usar)
survdiff(y ~ id_cat1, data = bancomg)
survdiff(y ~ cor1, data = bancomg)
survdiff(y ~ rendacat, data = bancomg)
survdiff(y ~ servico, data = bancomg)
survdiff(y ~ atraso2, data = bancomg)
survdiff(y ~ estadio, data = bancomg)

## Curvas de sobreviv�ncia pelo Kaplan-Meier estratificado

KMidade <- survfit(y ~ id_cat1, data = bancomg)
print(KMidade, print.rmean=TRUE)
summary(KMidade)
plot(KMidade, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:3, main = "Faixa et�ria", mark.time = F)
legend("bottomleft", c("<50","50-69",">=70"), lty = c(1:3), bty = "n")

KMcor <- survfit(y ~ cor1, data = bancomg)
print(KMcor, print.rmean=TRUE)
summary(KMcor)
plot(KMcor, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:2, main = "Ra�a/cor", mark.time = F)
legend("bottomleft", c("Branca", "N�o branca"), lty = c(1:2), bty = "n")

KMrenda <- survfit(y ~ rendacat, data = bancomg)
print(KMrenda, print.rmean=TRUE)
summary(KMrenda)
plot(KMrenda, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:4, main = "Renda", mark.time = F)
legend("bottomleft", c("Alta", "M�dia","Baixa","Muito baixa"), lty = c(1:4), bty = "n")

KMservico <- survfit(y ~ servico, data = bancomg)
print(KMservico, print.rmean=TRUE)
summary(KMservico)
plot(KMservico, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:2, main = "Setor", mark.time = F)
legend("bottomleft", c("Privado","P�blico"), lty = c(1:2), bty = "n")

KMatraso2 <- survfit(y ~ atraso2, data = bancomg)
print(KMatraso2, print.rmean=TRUE)
summary(KMatraso2)
plot(KMatraso2, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:2, main = "Atraso no tratamento", mark.time = F)
legend("bottomleft", c("N�o", "Sim"), lty = c(1:2), bty = "n")

KMestadio <- survfit(y ~ estadio, data = bancomg)
print(KMestadio, print.rmean=TRUE)
summary(KMestadio)
plot(KMestadio, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:4, main = "Estadiamento", mark.time = F)
legend("bottomleft", c("I","II","III","IV"), lty = c(1:4), bty = "n")

### MODELO SEMIPARAM�TRICO DE SOBREVIV�NCIA - RISCOS PROPORCIONAIS DE COX

library(survival)

## Criando objeto tipo sobreviv�ncia
y <- Surv(bancomg$tempo10,bancomg$status10)

## Modelos simples
summary(coxph(y ~ id_cat1, data = bancomg))
summary(coxph(y ~ cor1, data = bancomg))
summary(coxph(y ~ rendacat, data = bancomg))
summary(coxph(y ~ servico, data = bancomg))
summary(coxph(y ~ atraso2, data = bancomg))
summary(coxph(y ~ estadio, data = bancomg))

## Modelos m�ltiplos

# criando banco com dados completos para cor e renda
completo <- bancomg[complete.cases(bancomg[,c("cor1","rendacat")]),]

# modelos de Cox
cox1 <- coxph(y ~ cor1 + idade, data=completo)
summary(cox1)
cox.zph(cox1)

cox2 <- coxph(y ~ cor1 + idade + rendacat, data=completo)
summary(cox2)
cox.zph(cox2)
anova(cox1,cox2)

cox3 <- coxph(y ~ cor1 + idade + rendacat + servico, data=completo)
summary(cox3)
cox.zph(cox3)
anova(cox2,cox3)

cox4 <- coxph(y ~ cor1 + idade + rendacat + atraso2, data=completo)
summary(cox4)
cox.zph(cox4)
anova(cox2,cox4)

cox5 <- coxph(y ~ cor1 + idade + rendacat + estadio, data=completo)
summary(cox5)
cox.zph(cox5)
anova(cox2,cox5)

# testando intera��o entre as vari�veis

summary(coxph(y ~ cor1*idade, data=completo))
summary(coxph(y ~ cor1*rendacat, data=completo))
summary(coxph(y ~ cor1*jf, data=completo))
summary(coxph(y ~ cor1*servico, data=completo))
summary(coxph(y ~ cor1*atraso2, data=completo))
summary(coxph(y ~ cor1*estadio, data=completo))

summary(coxph(y ~ rendacat*idade, data=completo))
summary(coxph(y ~ rendacat*cor1, data=completo))
summary(coxph(y ~ rendacat*jf, data=completo))
summary(coxph(y ~ rendacat*servico, data=completo))
summary(coxph(y ~ rendacat*atraso2, data=completo))
summary(coxph(y ~ rendacat*estadio, data=completo))

## An�lise de res�duos do modelo cox5

# Calculando os res�duos de Schoenfeld (testa a correla��o linear global do modelo e de cada vari�vel)
zph5 <- cox.zph(cox5)
zph5

# Gr�fico dos Res�duos Schoenfeld de cada vari�vel do modelo cox5
par(mfrow=c(4,2),mar=c(2,2,2,2))
plot(zph5[1], main = "Cor")
abline(h = cox5$coef[1], lty = 3)
plot(zph5[2], main = "Idade")
abline(h = cox5$coef[2], lty = 3)
plot(zph5[3], main = "Renda m�dia")
abline(h = cox5$coef[3], lty = 3)
plot(zph5[4], main = "Renda baixa")
abline(h = cox5$coef[4], lty = 3)
plot(zph5[5], main = "Renda muito baixa")
abline(h = cox5$coef[5], lty = 3)
plot(zph5[6], main = "Estadio II")
abline(h = cox5$coef[6], lty = 3)
plot(zph5[7], main = "Estadio III")
abline(h = cox5$coef[7], lty = 3)
plot(zph5[8], main = "Estadio IV")
abline(h = cox5$coef[8], lty = 3)

# Avalia��o da exist�ncia de pontos aberrantes: res�duos deviance (fora do intervalo -2 a 2)
par(mfrow = c(1,1))
res.dev5 <- resid(cox5,type="deviance")
plot(res.dev5, xlab = "Indice", ylab = "Res�duo deviance", main= "Res�duos deviance")
abline(h=2,lty=3)
res.dev5[res.dev5>=2] # identifica indiv�duos com res�duos >=2

# Res�duo escore - pontos influentes 
res.esco5 <- resid(cox5,type="dfbetas")
plot(completo$cor1, res.esco5[,1], xlab = "Ra�a/cor", ylab = "Res�duos")
plot(completo$idade, res.esco5[,2], xlab = "Idade", ylab = "Res�duos")
plot(completo$rendacat, res.esco5[,3], xlab = "Renda m�dia", ylab = "Res�duos")
plot(completo$rendacat, res.esco5[,4], xlab = "Renda baixa", ylab = "Res�duos")
plot(completo$rendacat, res.esco5[,5], xlab = "Renda muito baixa", ylab = "Res�duos")
plot(completo$estadio, res.esco5[,6], xlab = "Estadio II", ylab = "Res�duos")
plot(completo$estadio, res.esco5[,7], xlab = "Estadio III", ylab = "Res�duos")
plot(completo$estadio, res.esco5[,8], xlab = "Estadio IV", ylab = "Res�duos")


### DECOMPOSI��O DE EFEITOS EM AN�LISE DE SOBREVIV�NCIA

# transformando algumas vari�veis
d = completo

d$cor2[d$cor1=="branca"] <- 0
d$cor2[d$cor1=="nao branca"] <- 1

d$estadio2[d$estadio=="I"] <- 1
d$estadio2[d$estadio=="II"] <- 2
d$estadio2[d$estadio=="III"] <- 3
d$estadio2[d$estadio=="IV"] <- 4

d$jf2[d$jf=="JF"] <- 0
d$jf2[d$jf=="Outras"] <- 1

# efeito total
library(survival)
TE = coxph(Surv(tempo10,status10) ~ cor2 + idade + rendacat + cluster(jf2), data=d)
summary(TE)

# efeito da exposi��o no mediador estadio
library(VGAM)
Mestadio = vglm(estadio2 ~ cor2+idade+rendacat,family=propodds,data=d)
summary(Mestadio)

# decomposi��o dos efeitos
doEffectDecomp = function(d)
{
  # Step 1: Replicate exposure variable, predict mediators
  d$cor2Temp = d$cor2
  Mestadio = vglm(estadio2 ~ cor2Temp+idade+rendacat,family=propodds,data=d)
  # Step 2: Replicate data with different exposures
  d1 = d2 = d
  d1$cor2star = 1
  d2$cor2star = 0
  newd = rbind(d1, d2)
  # Step 3: Compute weights for estadio
  newd$cor2Temp = newd$cor2
  direct = predict(Mestadio, newdata=newd, 
                   type='response')[cbind(1:nrow(newd),newd$estadio2)]
  newd$cor2Temp = newd$cor2star
  indirect = predict(Mestadio, newdata=newd, 
                   type='response')[cbind(1:nrow(newd),newd$estadio2)]
  newd$W = indirect/direct
  # Step 4: Weighted Cox Model
  cox = coxph(Surv(tempo10,status10) ~ cor2 + cor2star + idade +
                rendacat, weight=W, data=newd)
  # Return value: Estimates for total, direct, indirect effects
  TE = exp(sum(coef(cox)[c('cor2', 'cor2star')]))
  DE = exp(unname(coef(cox)['cor2']))
  IE = exp(sum(coef(cox)[c('cor2star')]))
  PM = log(IE) / log(TE)
  return(c(exp(coef(cox)), TE=TE, DE=DE, IE=IE, PM=PM))
}

doEffectDecomp(d)

# intervalos de confian�a por bootstrap
CSamp = function(d)
{
  s = sample(unique(d$jf2), replace=TRUE)
  return(do.call('rbind', lapply(s, function(x) d[d$jf2 == x, ])))
}
HRs = replicate(10000, doEffectDecomp(CSamp(d)))
apply(HRs, 1, quantile, c(0.025, 0.975))


### an�lise de sensibilidade

## intera��o
library(VGAM)
summary(vglm(estadio2 ~ cor2*(idade + rendacat),family=propodds, data=d))

## erros de classifica��o

# se todos os casos com ra�a/cor ignorada = branca
# (mudar o banco e rodar novamente os c�digos anteriores)
bancomg$cor2[bancomg$cor1=="branca"] <- 0
bancomg$cor2[is.na(bancomg$cor1)] <- 0
bancomg$cor2[bancomg$cor1=="nao branca"] <- 1

completo <- bancomg[complete.cases(bancomg[,c("cor2","rendacat")]),]
d = completo

d$estadio2[d$estadio=="I"] <- 1
d$estadio2[d$estadio=="II"] <- 2
d$estadio2[d$estadio=="III"] <- 3
d$estadio2[d$estadio=="IV"] <- 4

d$jf2[d$jf=="JF"] <- 0
d$jf2[d$jf=="Outras"] <- 1

# se todos os casos com ra�a/cor ignorada = n�o branca
bancomg$cor2[bancomg$cor1=="branca"] <- 0
bancomg$cor2[is.na(bancomg$cor1)] <- 1
bancomg$cor2[bancomg$cor1=="nao branca"] <- 1

completo <- bancomg[complete.cases(bancomg[,c("cor2","rendacat")]),]
d = completo

d$estadio2[d$estadio=="I"] <- 1
d$estadio2[d$estadio=="II"] <- 2
d$estadio2[d$estadio=="III"] <- 3
d$estadio2[d$estadio=="IV"] <- 4

d$jf2[d$jf=="JF"] <- 0
d$jf2[d$jf=="Outras"] <- 1

# se todos os casos com renda ignorada = alta
bancomg$rendacat[is.na(bancomg$rendacat)] <- "a_alta"

completo <- bancomg[complete.cases(bancomg[,c("cor1","rendacat")]),]
d = completo

d$cor2[d$cor1=="branca"] <- 0
d$cor2[d$cor1=="nao branca"] <- 1

d$estadio2[d$estadio=="I"] <- 1
d$estadio2[d$estadio=="II"] <- 2
d$estadio2[d$estadio=="III"] <- 3
d$estadio2[d$estadio=="IV"] <- 4

d$jf2[d$jf=="JF"] <- 0
d$jf2[d$jf=="Outras"] <- 1

# se todos os casos com renda ignorada = muito baixa
bancomg$rendacat[is.na(bancomg$rendacat)] <- "d_mbaixa"

completo <- bancomg[complete.cases(bancomg[,c("cor1","rendacat")]),]
d = completo

d$cor2[d$cor1=="branca"] <- 0
d$cor2[d$cor1=="nao branca"] <- 1

d$estadio2[d$estadio=="I"] <- 1
d$estadio2[d$estadio=="II"] <- 2
d$estadio2[d$estadio=="III"] <- 3
d$estadio2[d$estadio=="IV"] <- 4

d$jf2[d$jf=="JF"] <- 0
d$jf2[d$jf=="Outras"] <- 1


### Comparando duas curvas de sobreviv�ncia usando o tempo m�dio restrito

library(survival)
library(survRM2)

# transformando as vari�veis categ�ricas em num�ricas
completo$cor2[completo$cor1=="branca"] <- 0
completo$cor2[completo$cor1=="nao branca"] <- 1

completo$rendacat2[completo$rendacat=="a_alta"] <- 0
completo$rendacat2[completo$rendacat=="b_media"] <- 1
completo$rendacat2[completo$rendacat=="c_baixa"] <- 2
completo$rendacat2[completo$rendacat=="d_mbaixa"] <- 3

completo$jf2[completo$jf=="JF"] <- 0
completo$jf2[completo$jf=="Outras"] <- 1

completo$servico2[completo$servico=="privado"] <- 0
completo$servico2[completo$servico=="publico"] <- 1

completo$atraso3[completo$atraso2=="nao"] <- 0
completo$atraso3[completo$atraso2=="sim"] <- 1

completo$estadio2[completo$estadio=="I"] <- 0
completo$estadio2[completo$estadio=="II"] <- 1
completo$estadio2[completo$estadio=="III"] <- 2
completo$estadio2[completo$estadio=="IV"] <- 3

summary(completo)

# Estimando a curva de sobreviv�ncia pelo Kaplan-Meier
y <- Surv(completo$tempo10,completo$status10)
KM <- survfit(y ~ 1, data=completo)
print(KM, print.rmean=TRUE)

# Restricted mean survival time (RMST) and restricted mean time lost (RMTL)
# (�rea sob a curva da fun��o de sobrevida)
time   = completo$tempo10
status = completo$status10
arm    = completo$cor2
obj1 <- rmst2(time, status, arm)
obj1
plot(obj1, xlab="Years", ylab="Probability")

# fazendo a an�lise ajustada para idade (confundimento)
x2 = completo$idade
head(x2)
obj2 <- rmst2(time, status, arm, covariates=x2)
obj2

# fazendo a an�lise ajustada para idade e renda (confundimento)
x3 = completo[,c(2,29)]
head(x3)
obj3 <- rmst2(time, status, arm, covariates=x3)
obj3

# fazendo a an�lise ajustada para idade e renda (confundimento) e servi�o (media��o)
x4 = completo[,c(2,29,31)]
head(x4)
obj4 <- rmst2(time, status, arm, covariates=x4)
obj4

# fazendo a an�lise ajustada para idade e renda (confundimento) e atraso (media��o)
x5 = completo[,c(2,29,32)]
head(x5)
obj5 <- rmst2(time, status, arm, covariates=x5)
obj5

# fazendo a an�lise ajustada para idade e renda (confundimento) e estadiamento (media��o)
x6 = completo[,c(2,29,33)]
head(x6)
obj6 <- rmst2(time, status, arm, covariates=x6)
obj6

# comparando a diferen�a de tempo de sobreviv�ncia m�dio
# por ra�a/cor bruto e ajustado pelas covari�veis
obj1$unadjusted.result[1]
obj2$RMST.difference.adjusted[2]
obj3$RMST.difference.adjusted[2]
obj4$RMST.difference.adjusted[2]
obj5$RMST.difference.adjusted[2]
obj6$RMST.difference.adjusted[2]
