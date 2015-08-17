import graphics
from graphics import *

import random

# --------------------- Edge ----------------------
class Edge:
    v1 = 0
    v2 = 0
    w = 0
    shm = 0
    
    def __init__(self, v1, v2, w, shm):
        self.v1 = v1
        self.v2 = v2
        self.w = w
        self.shm = shm

# --------------------- Graph ---------------------
class Graph:
    edges = list()
    num_vertices = 0
    name = ""

    def __init(self):
        self.edges = list()

    def ParseNumVertices(self, first_line):
        self.num_vertices = int(first_line.strip()[2:])

    def ExtractFromFile(self, filename):
        self.name = filename
        fhandler = open(filename, "r")
        lines = fhandler.readlines()
        self.ParseNumVertices(lines[0])
        for i in range(1, len(lines)):
            splits = lines[i].strip().split()
            edge = Edge(int(splits[0]), int(splits[1]), int(splits[2]), int(splits[3]))
            self.edges.append(edge)
        fhandler.close()

# --------------------- Permutation ----------------
class Permutation:
    direct = list()
    reverse = list()

    def __init__(self, num_elems):
        self.direct = [0] * num_elems
        self.reverse = [0] * num_elems

    def ExtractFromFile(self, filename):
        fhandler = open(filename, "r")
        lines = fhandler.readlines()
        for i in range(0, len(lines)):
            num1 = i
            num2 = int(lines[i].strip())
            self.direct[num1] = num2
            self.reverse[num2] = num1
        fhandler.close()

# --------------------- Decomposition --------------
class Decomposition:
    vertex_color = list()
    class_color = dict()
    num_classes = 0

    def __init__(self):
        self.vertex_color = list()
        self.class_color = dict()
        self.num_classes = 0

    def ExtractFromFile(self, filename):
        fhandler = open(filename, 'r')
        lines = fhandler.readlines()
        for l in lines:
            set_id = int(l.strip())
            self.vertex_color.append(set_id)
            self.class_color[set_id] = 'red'
        self.num_classes = len(self.class_color)
        print "Decomposition into " + str(self.num_classes) + " classes were extracted from " + filename

    def AssignColors(self):
        index = 6
        for set_id in self.class_color:
            self.class_color[set_id] = dec_color_dict[index % len(dec_color_dict)]
            index += 1

# --------------------- drawing utils ---------------------
shm_color_dict = {0: color_rgb(41, 209, 37), 1: color_rgb(241, 143, 30), 2: color_rgb(239, 122, 26), 3: color_rgb(237, 102, 22), 4: color_rgb(235, 81, 19), 5: color_rgb(234, 61, 15), 6: color_rgb(232, 40, 12), 7: color_rgb(230, 20, 8), 8: color_rgb(229, 0, 5)}

weight_color_dict = {1: color_rgb(91, 255, 56), 2: color_rgb(65, 232, 46), 3: color_rgb(41, 209, 37), 4: color_rgb(29, 186, 37), 5: color_rgb(22, 163, 40), 6: color_rgb(16, 140, 42), 7: color_rgb(11, 117, 41), 8: color_rgb(7, 95, 39)}

dec_color_dict = {0: color_rgb(255, 0, 0), 1: color_rgb(255, 140, 0), 2: color_rgb(255, 255, 0), 3: color_rgb(0, 238, 118), 4: color_rgb(0, 191, 255), 5: color_rgb(0, 0, 255), 6: color_rgb(238, 0, 238)}

def DrawRectOnImg(img, x, y, x_size, y_size, color):
    for i in range(x, x + x_size):
        for j in range(y, y + y_size):
            img.setPixel(i, j, color)

def FillImage(img, img_size, color):
    for i in range(0, img_size):
        for j in range(0, img_size):
            img.setPixel(i, j, color)

def GetColorForValue(mode, edge, decomposition):
    if mode == 'weight':
        return weight_color_dict[edge.w]
    if mode == 'shm':
        return shm_color_dict[edge.shm]
    if decomposition.vertex_color[edge.v1] == decomposition.vertex_color[edge.v2]:
        return decomposition.class_color[decomposition.vertex_color[edge.v1]]
    return "grey"

def GetOutputFigName(graph, mode, max_tau, output_dir):
    output_file = os.path.join(output_dir, os.path.basename(graph.name) + "_tau_" + str(max_tau))
    if mode == 'weight':
        return output_file + "_weight.ppm"
    if mode == 'shm':
        return output_file + "_shm.ppm"
    return output_file + "_dec.ppm"

# --------------------- graph drawing ---------------------
def Draw_GraphSizeLessWindow(graph, perm, image, mode, max_tau, decomposition, window_size):
    rect_size = window_size / graph.num_vertices + 1
    for edge in graph.edges:
        if edge.w <= max_tau:
            x = int(float(perm.direct[edge.v1]) / graph.num_vertices * window_size)
            y = int(float(perm.direct[edge.v2]) / graph.num_vertices * window_size)
            DrawRectOnImg(image, x, y, rect_size, rect_size, GetColorForValue(mode, edge, decomposition))
            DrawRectOnImg(image, y, x, rect_size, rect_size, GetColorForValue(mode, edge, decomposition))
    if mode == 'dec':
        for i in range(0, graph.num_vertices):
            x = int(float(i) / graph.num_vertices * window_size)
            DrawRectOnImg(image, x, x, rect_size, rect_size, GetColorForValue(mode, Edge(perm.reverse[i], perm.reverse[i], 0, 0), decomposition))
    
def Draw_GraphSizeGreaterWindow(graph, perm, image, mode, max_tau, decomposition, window_size):
    for edge in graph.edges:
        if edge.w <= max_tau:
            x = int(float(perm.direct[edge.v1]) / graph.num_vertices * window_size)
            y = int(float(perm.direct[edge.v2]) / graph.num_vertices * window_size)
            image.setPixel(x, y, GetColorForValue(mode, edge, decomposition))
            image.setPixel(y, x, GetColorForValue(mode, edge, decomposition))

def DrawGraph(graph, perm, mode, max_tau, output_dir, decomposition, window_size = 600):
    image = Image(Point(0, 0), window_size, window_size)
    FillImage(image, window_size, "white")
    if graph.num_vertices <= window_size:
        Draw_GraphSizeLessWindow(graph, perm, image, mode, max_tau, decomposition, window_size)
    else:
        Draw_GraphSizeGreaterWindow(graph, perm, image, mode, max_tau, decomposition, window_size)
    image_fname = GetOutputFigName(graph, mode, max_tau, output_dir)
    print "Figure for mode " + mode + " and max_tau " + str(max_tau) + " was written to " + image_fname
    image.save(image_fname)

# --------------------- outline ---------------------
def DrawMatrices(graph_file, perm_file, mode, dec_file):
    graph = Graph()
    graph.ExtractFromFile(graph_file)

    perm = Permutation(graph.num_vertices)
    perm.ExtractFromFile(perm_file)

    decomposition = Decomposition()
    if dec_file != "":
        decomposition.ExtractFromFile(dec_file)
        decomposition.AssignColors()

    output_dir = os.path.basename(graph_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for tau in range(3, 9):
        DrawGraph(graph, perm, mode, tau, output_dir, decomposition)

    return output_dir

def main():
    graph_file = sys.argv[1]
    perm_file = sys.argv[2]
    mode = sys.argv[3] # weight or shm or dec
    dec_file = ""
    if mode == 'dec':
        dec_file = sys.argv[4]
    DrawMatrices(graph_file, perm_file, mode, dec_file)

if __name__ == "__main__":
    main()
