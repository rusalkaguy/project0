# map vs for

input_arr = [1,2,3,4]

print "# simple for loop for squares"
for x in input_arr:
	print x,x**2

print "# for loop with function"
def p(x):
	print x, x**2
	return str(x**2)

for x in input_arr:
	p(x)

print "# map with function"
result_arr = map(p,input_arr)
print "result_arr=",result_arr
print "\n".join(result_arr)

print "# map with anonymous function"
result_arr = map(lambda x: str(x**2), input_arr)
print "result_arr=",result_arr
print "\n".join(result_arr)

