#################################################################################
#   Arquivo: Gaussiana-Inversa.R
#  
#   USO: Simulacao de Monte Carlo dos estimadores de maxima
#   verossimilhanca dos parametros da distribuiCao Gaussiana Inversa, IG(m,l).
#   
#   OTIMIZACAO NAO-LINEAR: Feita por com gradiente analitico, método BFGS
#
#   ARGUMENTOS DA FUNCAO: nobs -> numero de observacoes testados foram 
#   (T =25,50,75,100,250,400);
#                         nrep -> numero de replicas de Monte Carlo (R = 10000);
#                         semente -> semente do gerador (2000);
#                         shape1 -> valor de mu
#                         shape2 -> valor de lambda
#
#   Autora: Gabrielle Carvalho
#
#   Data: 07/07/2022
#
############################################################################

# pacote para geracao aleatoria
install.packages(	"statmod")
library("statmod")

# Criação da funcao
gaussInv.mc.bfgs = function(nrep=10000, nobs=100, semente=2000,
                          shape1 = 1.5, shape2= 2)
{  
  # função de log-verossimilhança
  logLikGaussInv = function(theta){
    a = abs(theta[1]) 
    b = abs(theta[2])
    logLik =  sum( 0.5*log(b) - 1.5*log(y) - (b*(y - a)^2 )/(2* a^2 *y) )
    return(logLik)
  }


  # funcao escore (gradiente)
  scoreFn = function(theta){
    a = abs(theta[1]);  b = abs(theta[2])
    cbind( sum((b/(a^3*y)) * (a*(y-a) + (y-a)^2 )),
            sum(1/(2*b)  - ((y-a)^2 )/(2* a^2 *y) ))
  }

  # início da cronometragem
  tempo.inicio = Sys.time()

  # vetores para armanezar as estimativas
  emva = rep(0, nrep)
  emvb = rep(0, nrep)

  set.seed(2000) # semente do gerador
  contadorFalhas = 0 # contador de falhas

  # laço de Monte Carlo
  i = 1
  while(i <= nrep){
  
    # amostra gerada
    y = rinvgauss(nobs, mean = shape1, dispersion = 1/shape2)
  
  
    # chute inicial
    shoot1 = cbind(1.3)
    shoot2 = cbind(1.7)
  
    # maximização da log-verossimilhança
    ir = optim(c(shoot1, shoot2), logLikGaussInv, gr=scoreFn, method="BFGS",
               control=list(fnscale=-1, reltol = 1e-12))
  
    # checagem de convergencia
    if(ir$convergence == 0){
      emva[i] = ir$par[1]
      emvb[i] = ir$par[2]
      i = i + 1
    }
    else{
      contadorFalhas = contadorFalhas + 1
    }
  } # fim do laco de Monte Carlo

  # estimativas
  amedio = mean(emva)
  bmedio = mean(emvb)
  avies = amedio - shape1
  bvies = bmedio - shape2
  aviesrel = avies/shape1*100
  bviesrel = bvies/shape2*100
  amax = max(emva)
  bmax = max(emvb)
  amin = min(emva)
  bmin = min(emvb)
  
  # erros quadraticos madios
  aeqm = avies^2 + var(emva)
  beqm = bvies^2 + var(emvb)
  
  
  # resultados armazenados em matriz
  mResultados = matrix(c(amedio, bmedio, avies,
                         bvies, aviesrel, bviesrel,
                         aeqm, beqm, amax, bmax,
                         amin,bmin), 6, 2, byrow=TRUE)
  rownames(mResultados) = c("média", "viés",
                            "viés rel.", "eqm", "max", "min")
  colnames(mResultados) = c("evm a", "emv b")
  
  # calculo do tempo de execucao
  tempo.fim = Sys.time()
  tempo.exec = tempo.fim - tempo.inicio
  
  # resultados colecionados em uma lista
  resultado = list(nobs=nobs, nrep=nrep, semente=semente, a=shape1,
                   b=shape2, falhas=contadorFalhas, resultados=mResultados,
                   horario=tempo.inicio, tempoexec=tempo.exec)
  return(resultado)
}
#res50 = gaussInv.mc.bfgs(nobs=50)
#res50$resultados
sessioninfo() 