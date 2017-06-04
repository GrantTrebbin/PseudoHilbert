import PseudoHilbert
#import time
#import cProfile
import svgwrite


class Diagram:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.dwg = svgwrite.Drawing(profile='full',
                                    size=(str(width) + 'px',
                                          str(height) + 'px'),
                                    viewBox=('0 0 ' + str(width) +
                                             ' ' + str(height)))
        self.dwg.add(self.dwg.rect(insert=(0, 0),
                                   size=(width, height),
                                   fill='rgb(255,255,255)'))

rectangle_width = 23
rectangle_height = 17

diagram = Diagram(rectangle_width, rectangle_height)

PsH = PseudoHilbert.PseudoHilbert(rectangle_width, rectangle_height)
path = PsH.index_to_coordinate
coord_adjusted_path = [[coord[0] + 0.5, rectangle_height - coord[1] - 0.5]
                       for coord in path]

diagram.dwg.add(diagram.dwg.polyline(coord_adjusted_path,
                                     stroke="rgb(0,0,0)",
                                     stroke_linejoin='round',
                                     stroke_width=0.3,
                                     fill='none'))

diagram.dwg.filename = format(0, '04') + '.svg'
diagram.dwg.save()




#print(time.time())
#cProfile.run(
#    "test = PseudoHilbert.PseudoHilbert(191, 192)")
#cProfile.run(
#    "test = PseudoHilbert(11*270, 11*190)")
#print(time.time())
