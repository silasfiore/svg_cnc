from subprocess import call
from xml.dom import minidom


# call dat2svg, name the two curves left0, left1, left2 ... and right... or sth to identify them.
call(["./build/dat2svg","miley.dat","-c","120","-y","30","-x","130","-a","180","-r"])
call(["./build/dat2svg","sd7055.dat","-c","120","-y","30","-x","130","-a","178","-r"])




# add entry and exit curves to trajectory ( this can easily be done here!)


# read the SVG file
doc = minidom.parse('miley.svg')
path1 = [path.getAttribute('d') for path in doc.getElementsByTagName('path')]
doc.unlink()

doc = minidom.parse('sd7055.svg')
path2 = [path.getAttribute('d') for path in doc.getElementsByTagName('path')]
doc.unlink()


path1.insert(0,"L 0,30")
path2.insert(0,"L 0,30")


# print the line draw commands
#for path_string in path_strings:



path1[1] = "M 0,30 " + path1[1].replace("M","L",1)
path2[1] = "M 0,30 " + path2[1].replace("M","L",1)

call(["./build/prog", path1[0]+path1[1]+path1[2],path2[0]+path2[1]+path2[2]])

for seg1,seg2 in zip(path1,path2):
    call(["./build/prog",seg1,seg2])


#call(["./build/prog", T3, T4 ])