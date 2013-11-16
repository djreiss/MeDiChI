library( MeDiChI )
data( "halo.lowres", package="MeDiChI" )
cm.fit <- chip.deconv( data.halo.lowres, where="Chr", fit.res=10, center=650000, wind=20000,
                      max.steps=100, n.boot=10, kernel=kernel.halo.lowres, verbose=TRUE, boot.sample.opt="case" )
plot( cm.fit, plot.genes=TRUE )
