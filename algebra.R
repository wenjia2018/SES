foldl = function(x, ...) reduce(x, ~str_c("f(", .x, " , ", .y, ")"), .dir = "forward", ...)
foldr = function(x, ...) reduce(x, ~str_c("f(", .x, " , ", .y, ")"), .dir = "backward", ...)

scanl = function(x, ...) accumulate(x, ~str_c("f(", .x, " , ", .y, ")"), .dir = "forward", ...)
scanr = function(x, ...) accumulate(x, ~str_c("f(", .x, " , ", .y, ")"), .dir = "backward", ...)

# recurse into the second (right) argument of f
foldr(1:10, .init = 0)
scanr(1:10, .init = 0)

# here .init is the right identity (the second identity), and acts as 

# recurse into the first argument of f
foldl(1:10, .init = 0)
scanl(1:10, .init = 0)
# here .init is the left identity


# infix operator notation for the function
foldr2 = function(x, ...) reduce(x, ~str_c("(", .x, " f ", .y, ")"), .dir = "backward", ...)
foldr2(1:10)

# Like a right fold, a left fold cannot perform magic and go to the end of the
# list instantly; it must start from the beginning of the list. However, the
# parentheses dictate how our code evaluates. The type of the argument to the
# folding function changes in addition to the associativity:

# The fundamental way to think about evaluation in Haskell is as substitution.
# When we use a right fold on a list with the function 𝑓 and start value 𝑧,
# we’re, in a sense, replacing the cons constructors - in R, these are the
# commas - with our folding function and the empty list constructor with our
# start value z

# a la Haskell
cons = function(x) reduce(x, function(x, y) c(x, y), .init = list()) 

cons(c(1, 3, 4))


# integers under + form a group, hence a monoid
reduce(1:10, `+`)
accumulate(1:10, `+`, .dir = "backward")
sum(1:10)


# permutation group examples, but why does this appear commutative?
my_perms = replicate(10, sample(5), simplify = FALSE)

reduce(my_perms, PerMallows::compose, .dir = "backward")
reduce(my_perms, PerMallows::compose, .dir = "forward")

accumulate(my_perms, PerMallows::compose, .dir = "backward")
accumulate(my_perms, PerMallows::compose, .dir = "forward")

#

rduce(list(mean, sqrt), compose, .dir = "forward")(1:10)
mean(sqrt(1:10))

