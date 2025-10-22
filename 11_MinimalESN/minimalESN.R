# A minimalistic Echo State Networks demo with Mackey-Glass (delay 17) data 
# in "plain" R.
# by Mantas Lukosevicius 2012-2018
# http://mantas.info

# load the data
trainLen = 2000
testLen = 2000
initLen = 100
data = as.matrix(read.table('MackeyGlass_t17.txt'))

# plot some of it
while( dev.cur() != 1 ) dev.off() # close all previous plots
dev.new()
plot(data[1:1000],type='l')
title(main='A sample of data')

# generate the ESN reservoir
inSize = outSize = 1
resSize = 1000
a = 0.3 # leaking rate

set.seed(42)
Win = matrix(runif(resSize*(1+inSize),-0.5,0.5),resSize)
W = matrix(runif(resSize*resSize,-0.5,0.5),resSize)

# normalizing and setting spectral radius
cat('Computing spectral radius...')
rhoW = abs(eigen(W,only.values=TRUE)$values[1])
print('done.')
W = W * 1.25 / rhoW

# allocated memory for the design (collected states) matrix
X = matrix(0,1+inSize+resSize,trainLen-initLen)
# set the corresponding target matrix directly
Yt = matrix(data[(initLen+2):(trainLen+1)],1)

# run the reservoir with the data and collect X
x = rep(0,resSize)
for (t in 1:trainLen){
    u = data[t]
    x = (1-a)*x + a*tanh( Win %*% rbind(1,u) + W %*% x )
    if (t > initLen)
        X[,t-initLen] = rbind(1,u,x)
}

# train the output
reg = 1e-8  # regularization coefficient
X_T = t(X)
Wout = Yt %*% X_T %*% solve( X %*% X_T + reg*diag(1+inSize+resSize) )

# run the trained ESN in a generative mode. no need to initialize here, 
# because x is initialized with training data and we continue from there.
Y = matrix(0,outSize,testLen)
u = data[trainLen+1]
for (t in 1:testLen){
    x = (1-a)*x + a*tanh( Win %*% rbind(1,u) + W %*% x )
    y = Wout %*% rbind(1,u,x)
    Y[,t] = y
    # generative mode:
    u = y
    # this would be a predictive mode:
    #u = data[trainLen+t+1] 
}

# compute MSE for the first errorLen time steps
errorLen = 500
mse = ( sum( (data[(trainLen+2):(trainLen+errorLen+1)] - Y[1,1:errorLen])^2 )
	/ errorLen )
print( paste( 'MSE = ', mse ) )

# plot some signals
dev.new() 
plot( data[(trainLen+1):(trainLen+testLen+1)], type='l', col='green' )
lines( c(Y), col='blue' )
title(main=expression(paste('Target and generated signals ', bold(y)(italic(n)), 
	' starting at ', italic(n)==0 )))
legend('bottomleft',legend=c('Target signal', 'Free-running predicted signal'), 
	col=c('green','blue'), lty=1, bty='n' )

dev.new()
matplot( t(X[(1:20),(1:200)]), type='l' )
title(main=expression(paste('Some reservoir activations ', bold(x)(italic(n)))))

dev.new()
barplot( Wout )
title(main=expression(paste('Output weights ', bold(W)^{out})))