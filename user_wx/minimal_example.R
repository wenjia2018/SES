myfun = function(x,y,z){
  # out = NULL
  out$detail = x+y
  out$p= x+z
  out$other$mediation = max(z,out$detail)
  out
}

myfun(1,3,5)
