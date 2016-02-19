# learning regex by example
import re
   
#pat= '[^aeiouy]+([aeiouy]+)'  # grouping operator on vowels
pat= '([^aeiouy]+)[aeiouy]+'  # grouping operator on non-vowels
#pat= '[^aeiouy]+[aeiouy]+'  # no grouping operator on vowels
target = 'blair writes programs'

target='complement(join(48595..49823,49927..50128))'
pat='complement\((.*)\)'
print "target=", target
p = re.compile(pat)
#------------------------------------------
m = p.match( target )
print "pat=",pat, " match=", m
if m:
    print 'Match found. Span=', m.span(), ' Group(0)=', m.group(), 'Group(1)=', m.group(1)
else:
    print 'No match'
#------------------------------------------
s = p.search( target )
print "pat=",pat, " search=", s
if s:
    print 'search found. Span=', s.span(), ' Group=', s.group()
else:
    print 'No search hit'
#------------------------------------------
fa = p.findall( target )
print "pat=",pat, " findall=", fa
#------------------------------------------
hits = p.finditer(target)
for match in hits:
    print 'finditer found. Span=', match.span(), ' Group=', match.group()
#    print 'No search hit'