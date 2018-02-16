library(gramEvol)
library(xtable)
library(pracma)
library(stringr)
library(stringi)
#-------------------------------------Data---------------------------------------
#categories
category<-LETTERS[1:10]
#distribution of the error
error<-runif(10)
#Correlation between Y and error
y<-0.5*cos(error)

#------------------Design of Rules and BNF Grammar-----------------------------
system.rules.1<- list(expr = grule(op(expr, expr), 
                                   func(expr), 
                                   var),
                func = grule(sin, cos, log, sqrt), 
                op = grule(`+`, `-`, `*`),
                var = grule(error, 
                            n*error, 
                            n), 
                n = grule(0.25, 0.5, 0.75, 1))

grammar.1<- CreateGrammar(system.rules.1)
test.grammar<-summary(grammar.1)

system.rules.2 <- list(expr = gsrule("<var><op><var>"),
                op   = gsrule("+", "-", "*"),
                var  = gvrule(category))

grammar.2<- CreateGrammar(system.rules.2)
test.grammar.2<-summary(grammar.2)

#-------------------Evaluate the Design of the Grammar ------------------------

grammar.sequence <- c(0, 0, 1, 2)
grammar.expression<- GrammarMap(grammar.sequence, grammar.1, verbose = TRUE)
grammar.expressions <- GrammarRandomExpression(grammar.1, 5)

grammar.evaluations.df<-as.data.frame(EvalExpressions(grammar.expressions, data.frame(error=error)))

Design.grammar.df<-data.frame()
Design.grammar.df<-rbind(c("Grammar 1",
                           test.grammar$No.of.Unique.Expressions,
                           test.grammar$Tree.Depth),
                           test.grammar$Maximum.Sequence.Length)
colnames(Design.grammar.df)<-c("Grammar","Unique Expressions","Tree Depth","Max Sequence Length")

#---------- Grammatical Evolution and Other Search Algorithms--------------------------------------------

cost.k<-0.021
ge.solution.1 <- GrammaticalEvolution(grammar.1, 
                                      objective.functional,
                                      max.depth=GrammarGetDepth(grammar.1),
                                      seqLen=GrammarMaxSequenceLen(grammar.1,4,test.grammar$Start.Symbol),
                                      terminationCost = cost.k,
                                      popSize="auto",
                                      optimizer="auto",
                                      newPerGen="auto",
                                      elitism=2,
                                      iterations=500) 

ges.solution.1 <- GrammaticalExhaustiveSearch(grammar.1, 
                                              objective.functional, 
                                              max.depth=2,
                                              iterations=100,
                                              monitorFunc=algorithm.visualization,
                                              terminationCost = cost.k)

grs.solution.1 <- GrammaticalRandomSearch(grammar.1, 
                                          objective.functional,
                                          max.depth=2,
                                          iterations=100,
                                          terminationCost = cost.k,
                                          monitorFunc=algorithm.visualization)


ge.solution.1.best.expression <- ge.solution.1$best$expression
ges.solution.1.best.expression <- ges.solution.1$best$expression
grs.solution.1.best.expression <- grs.solution.1$best$expression

ge.test.df<-data.frame(category,
                       error, 
                       y, 
                       Test = 0.5*cos(error), 
                       GE = eval(ge.solution.1.best.expression))

ges.test.df<-data.frame(category,
                        error, 
                       y, 
                       Test = 0.5*error, 
                       GE = eval(ges.solution.1.best.expression))

grs.test.df<-data.frame(category,
                        error, 
                       y, 
                       Test = 0.5*error, 
                       GE = eval(grs.solution.1.best.expression))

ge.solution.expressions.df<-data.frame()
ge.solution.expressions.df<-rbind(c(ge.solution.1$settings$iterations,
                                    "GE",
                                    as.character(ge.solution.1.best.expression))
                                  )
colnames(ge.solution.expressions.df)<-c("Iterations","Algorithm","Best Expression")


#-----------------------------------Tables---------------------------------------

Table.1<-xtable(Design.grammar.df)
Table.2<-xtable(ge.test.df)
Table.3<-xtable(ges.test.df)
Table.4<-xtable(grs.test.df)
Table.5<-xtable(grammar.evaluations.df)

#-----------------------------------Figures-------------------------------------

Figure.1<-plot(seq(1:length(ge.test.df$Test)),ge.test.df$Test)
lines(ge.test.df$GE)

Figure.2<-hist(error)

#-----------------------------------References----------------------------------

Reference.1<-c("",
               "",
               "")

#-----------------------------------Function Library-----------------------------

grammar.check<-function(x,y,n)
{
  string <- NULL
  for(i in 1:n){
    string <- GrammarGetNextSequence(x,y)
    if (is.GrammarOverflow(string)) {
      break
    }
    expr <- GrammarMap(string, x)
    cat(string, " -> ", as.character(expr), "\n")
  }
}

test.grammar.check<-grammar.check(grammar.1,c(0,1),3)

objective.functional<-function(expression)
{
  result <- eval(expression) 
  if (any(is.nan(result))) 
    return(Inf) 
   res<-mean(log(1 + abs(y - result)))
   return(res)
}


algorithm.visualization <- function(results)
  { 
  #print(results) 
  } 

