
# -*- coding: utf-8 -*-


class Rectangle:
    def __init__(self):
        self.width = 0
        self.height = 0

    def setSize(self,size):
        self.width,self.height = size

    def getSize(self):
        return self.width, self.height

    size = property(getSize,setSize)

r = Rectangle()
r.width = 10
r.height = 5
print r.size
print r.size
