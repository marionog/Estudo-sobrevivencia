### ESTUDO DE SOBREVIVÊNCIA EM 10 ANOS DE
### MULHERES COM CÂNCER DE MAMA EM JUIZ DE FORA
### DIAGNOSTICADAS ENTRE 2003-2005

### ANALISES EXPLORATÓRIAS

## RESUMO NUMÉRICO DO BANCO
summary(bancomg)

## ANÁLISE COM FOCO NA VARIÁVEL RAÇA/COR (BRANCA/NÃO BRANCA)

# Relação entre cor e faixa etária
cor.fxidade <- table(bancomg$id_cat1,bancomg$cor1)
cor.fxidade
round(prop.table(cor.fxidade,1),3)
round(prop.table(cor.fxidade,2),3)
summary(cor.fxidade)
fisher.test(cor.fxidade)

# Relação entre cor e renda
cor.renda <- table(bancomg$rendacat,bancomg$cor1)
cor.renda
round(prop.table(cor.renda,1),3)
round(prop.table(cor.renda,2),3)
summary(cor.renda)
fisher.test(cor.renda)

# Relação entre cor e natureza do serviço
cor.servico <- table(bancomg$servico,bancomg$cor1)
cor.servico
round(prop.table(cor.servico,1),3)
round(prop.table(cor.servico,2),3)
summary(cor.servico)
fisher.test(cor.servico)

# Relação entre cor e atraso no início do tratamento
cor.atraso <- table(bancomg$atraso2,bancomg$cor1)
cor.atraso
round(prop.table(cor.atraso,1),3)
round(prop.table(cor.atraso,2),3)
summary(cor.atraso)
fisher.test(cor.atraso)

# Relação entre cor e estadiamento do tumor
cor.estadio <- table(bancomg$estadio,bancomg$cor1)
cor.estadio
round(prop.table(cor.estadio,1),3)
round(prop.table(cor.estadio,2),3)
summary(cor.estadio)
fisher.test(cor.estadio)

# Relação entre cor e óbito - 10 anos
cor.status10 <- table(bancomg$status10,bancomg$cor1)
cor.status10
round(prop.table(cor.status10,1),3)
round(prop.table(cor.status10,2),3)
summary(cor.status10)
fisher.test(cor.status10)


### ANÁLISE NÃO PARAMÉTRICA - KAPLAN-MEIER

library(survival)

# Criando objeto tipo sobrevivência
y <- Surv(bancomg$tempo10,bancomg$status10)

# Estimando a curva de sobrevivência pelo Kaplan-Meier
KM <- survfit(y ~ 1, data=bancomg)
print(KM, print.rmean=TRUE)
options(max.print=10000)
summary(KM)
plot(KM, xlab = "Tempo (dias)", ylab = "S(t)", main = "Geral", mark.time = F)

## Teste Log-rank
survdiff(y ~ id_cat1, data = bancomg)
survdiff(y ~ cor1, data = bancomg)
survdiff(y ~ rendacat, data = bancomg)
survdiff(y ~ servico, data = bancomg)
survdiff(y ~ atraso2, data = bancomg)
survdiff(y ~ estadio, data = bancomg)

## Curvas de sobrevivência pelo Kaplan-Meier estratificado

KMidade <- survfit(y ~ id_cat1, data = bancomg)
print(KMidade, print.rmean=TRUE)
summary(KMidade)
plot(KMidade, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:3, main = "Faixa etária", mark.time = F)
legend("bottomleft", c("<50","50-69",">=70"), lty = c(1:3), bty = "n")

KMcor <- survfit(y ~ cor1, data = bancomg)
print(KMcor, print.rmean=TRUE)
summary(KMcor)
plot(KMcor, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:2, main = "Raça/cor", mark.time = F)
legend("bottomleft", c("Branca", "Não branca"), lty = c(1:2), bty = "n")

KMrenda <- survfit(y ~ rendacat, data = bancomg)
print(KMrenda, print.rmean=TRUE)
summary(KMrenda)
plot(KMrenda, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:4, main = "Renda", mark.time = F)
legend("bottomleft", c("Alta", "Média","Baixa","Muito baixa"), lty = c(1:4), bty = "n")

KMservico <- survfit(y ~ servico, data = bancomg)
print(KMservico, print.rmean=TRUE)
summary(KMservico)
plot(KMservico, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:2, main = "Setor", mark.time = F)
legend("bottomleft", c("Privado","Público"), lty = c(1:2), bty = "n")

KMatraso2 <- survfit(y ~ atraso2, data = bancomg)
print(KMatraso2, print.rmean=TRUE)
summary(KMatraso2)
plot(KMatraso2, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:2, main = "Atraso no tratamento", mark.time = F)
legend("bottomleft", c("Não", "Sim"), lty = c(1:2), bty = "n")

KMestadio <- survfit(y ~ estadio, data = bancomg)
print(KMestadio, print.rmean=TRUE)
summary(KMestadio)
plot(KMestadio, conf.int = F, xlab = "Tempo (dias)", ylab = "S(t)", lty = 1:4, main = "Estadiamento", mark.time = F)
legend("bottomleft", c("I","II","III","IV"), lty = c(1:4), bty = "n")

### MODELO SEMIPARAMÉTRICO DE SOBREVIVÊNCIA - RISCOS PROPORCIONAIS DE COX

library(survival)

## Criando objeto tipo sobrevivência
y <- Surv(bancomg$tempo10,bancomg$status10)

## Modelos simples
summary(coxph(y ~ id_cat1, data = bancomg))
summary(coxph(y ~ cor1, data = bancomg))
summary(coxph(y ~ rendacat, data = bancomg))
summary(coxph(y ~ servico, data = bancomg))
summary(coxph(y ~ atraso2, data = bancomg))
summary(coxph(y ~ estadio, data = bancomg))

## Modelos múltiplos

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

# testando interação entre as variáveis

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

## Análise de resíduos do modelo cox5

# Calculando os resíduos de Schoenfeld (testa a correlação linear global do modelo e de cada variável)
zph5 <- cox.zph(cox5)
zph5

# Gráfico dos Resíduos Schoenfeld de cada variável do modelo cox5
par(mfrow=c(4,2),mar=c(2,2,2,2))
plot(zph5[1], main = "Cor")
abline(h = cox5$coef[1], lty = 3)
plot(zph5[2], main = "Idade")
abline(h = cox5$coef[2], lty = 3)
plot(zph5[3], main = "Renda média")
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

# Avaliação da existência de pontos aberrantes: resíduos deviance (fora do intervalo -2 a 2)
par(mfrow = c(1,1))
res.dev5 <- resid(cox5,type="deviance")
plot(res.dev5, xlab = "Indice", ylab = "Resíduo deviance", main= "Resíduos deviance")
abline(h=2,lty=3)
res.dev5[res.dev5>=2] # identifica indivíduos com resíduos >=2

# Resíduo escore - pontos influentes 
res.esco5 <- resid(cox5,type="dfbetas")
plot(completo$cor1, res.esco5[,1], xlab = "Raça/cor", ylab = "Resíduos")
plot(completo$idade, res.esco5[,2], xlab = "Idade", ylab = "Resíduos")
plot(completo$rendacat, res.esco5[,3], xlab = "Renda média", ylab = "Resíduos")
plot(completo$rendacat, res.esco5[,4], xlab = "Renda baixa", ylab = "Resíduos")
plot(completo$rendacat, res.esco5[,5], xlab = "Renda muito baixa", ylab = "Resíduos")
plot(completo$estadio, res.esco5[,6], xlab = "Estadio II", ylab = "Resíduos")
plot(completo$estadio, res.esco5[,7], xlab = "Estadio III", ylab = "Resíduos")
plot(completo$estadio, res.esco5[,8], xlab = "Estadio IV", ylab = "Resíduos")


### DECOMPOSIÇÃO DE EFEITOS EM ANÁLISE DE SOBREVIVÊNCIA

# transformando algumas variáveis
d = completo

d$cor2[d$cor1=="branca"] <- 0
d$cor2[d$cor1=="nao branca"] <- 1

d$estadio2[d$estadio=="I"] <- 0
d$estadio2[d$estadio=="II"] <- 0
d$estadio2[d$estadio=="III"] <- 1
d$estadio2[d$estadio=="IV"] <- 1

d$jf2[d$jf=="JF"] <- 0
d$jf2[d$jf=="Outras"] <- 1

# pacotes necessários
library(survival)

# decomposição dos efeitos
doEffectDecomp = function(d)
{
  # Step 1: Replicate exposure variable, predict mediators
  d$cor2Temp = d$cor2
  Mestadio = glm(estadio2 ~ cor2Temp+idade+rendacat,family=binomial,data=d)
  # Step 2: Replicate data with different exposures
  d1 = d2 = d
  d1$cor2star = d1$cor2
  d2$cor2star = !di$cor2
  newd = rbind(d1, d2)
  # Step 3: Compute weights for estadio
  newd$cor2Temp = newd$cor2
  w = predict(Mestadio, newdata=newd, type='response')
  direct = ifelse(newd$estadio2, w, 1-w)
  newd$cor2Temp = newd$cor2star
  w = predict(Mestadio, newdata=newd, type='response')
  indirect = ifelse(newd$estadio2, w, 1-w)
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

# intervalos de confiança por bootstrap
CSamp = function(d)
{
  s = sample(unique(d$jf2), replace=TRUE)
  return(do.call('rbind', lapply(s, function(x) d[d$jf2 == x, ])))
}
HRs = replicate(10000, doEffectDecomp(CSamp(d)))
apply(HRs, 1, quantile, c(0.025, 0.975))


### análise de sensibilidade

## interação entre variável de exposição e mediadora
summary(coxph(Surv(tempo10,status10) ~ cor2*estadio2 + idade + rendacat + cluster(jf2), data=d))

## erros de classificação

# se todos os casos com raça/cor ignorada = branca
# (mudar o banco e rodar novamente os códigos anteriores)
bancomg$cor2[bancomg$cor1=="branca"] <- 0
bancomg$cor2[is.na(bancomg$cor1)] <- 0
bancomg$cor2[bancomg$cor1=="nao branca"] <- 1
completo <- bancomg[complete.cases(bancomg[,c("cor2","rendacat")]),]

# se todos os casos com raça/cor ignorada = não branca
# (mudar o banco e rodar novamente os códigos anteriores)
bancomg$cor2[bancomg$cor1=="branca"] <- 0
bancomg$cor2[is.na(bancomg$cor1)] <- 1
bancomg$cor2[bancomg$cor1=="nao branca"] <- 1
completo <- bancomg[complete.cases(bancomg[,c("cor2","rendacat")]),]

# se todos os casos com renda ignorada = alta
# (mudar o banco e rodar novamente os códigos anteriores)
bancomg$rendacat[is.na(bancomg$rendacat)] <- "a_alta"
completo <- bancomg[complete.cases(bancomg[,c("cor1","rendacat")]),]

# se todos os casos com renda ignorada = muito baixa
# (mudar o banco e rodar novamente os códigos anteriores)
bancomg$rendacat[is.na(bancomg$rendacat)] <- "d_mbaixa"
completo <- bancomg[complete.cases(bancomg[,c("cor1","rendacat")]),]


### Comparando duas curvas de sobrevivência usando o tempo médio restrito

library(survival)
library(survRM2)

# transformando as variáveis categóricas em numéricas
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

# Estimando a curva de sobrevivência pelo Kaplan-Meier
y <- Surv(completo$tempo10,completo$status10)
KM <- survfit(y ~ 1, data=completo)
print(KM, print.rmean=TRUE)

# Restricted mean survival time (RMST) and restricted mean time lost (RMTL)
# (área sob a curva da função de sobrevida)
time   = completo$tempo10
status = completo$status10
arm    = completo$cor2
obj1 <- rmst2(time, status, arm)
obj1
plot(obj1, xlab="Years", ylab="Probability")

# fazendo a análise ajustada para idade (confundimento)
x2 = completo$idade
head(x2)
obj2 <- rmst2(time, status, arm, covariates=x2)
obj2

# fazendo a análise ajustada para idade e renda (confundimento)
x3 = completo[,c(2,29)]
head(x3)
obj3 <- rmst2(time, status, arm, covariates=x3)
obj3

# fazendo a análise ajustada para idade e renda (confundimento) e serviço (mediação)
x4 = completo[,c(2,29,31)]
head(x4)
obj4 <- rmst2(time, status, arm, covariates=x4)
obj4

# fazendo a análise ajustada para idade e renda (confundimento) e atraso (mediação)
x5 = completo[,c(2,29,32)]
head(x5)
obj5 <- rmst2(time, status, arm, covariates=x5)
obj5

# fazendo a análise ajustada para idade e renda (confundimento) e estadiamento (mediação)
x6 = completo[,c(2,29,33)]
head(x6)
obj6 <- rmst2(time, status, arm, covariates=x6)
obj6

# comparando a diferença de tempo de sobrevivência médio
# por raça/cor bruto e ajustado pelas covariáveis
obj1$unadjusted.result[1]
obj2$RMST.difference.adjusted[2]
obj3$RMST.difference.adjusted[2]
obj4$RMST.difference.adjusted[2]
obj5$RMST.difference.adjusted[2]
obj6$RMST.difference.adjusted[2]

