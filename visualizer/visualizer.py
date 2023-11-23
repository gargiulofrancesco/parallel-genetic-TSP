from graphics import *
import json


def draw_city(x,y,win):
    pt = Point(x, y)
    outcir = Circle(pt, 10)
    cir = Circle(pt, 4)
    cir.setFill('black')
    outcir.draw(win)
    cir.draw(win)


def draw_line(x1,y1,x2,y2,win):
    pt1 = Point(x1, y1)
    pt2 = Point(x2, y2)
    ln = Line(pt1, pt2)
    ln.setOutline('red')
    ln.setWidth(3)
    ln.draw(win)


def draw_text(curr_gen, curr_distance, win):
    x = 175
    y = 500
    gen = Text(Point(x, y), "Generation: " + str(curr_gen))
    gen.setSize(16)
    gen.setFace('courier')
    gen.setTextColor('black')
    gen.draw(win)

    dis = Text(Point(x+40,y+25), "Distance: " + str("%.2f" % curr_distance))
    dis.setSize(16)
    dis.setFace('courier')
    dis.setTextColor('black')
    dis.draw(win)


def clean_window(win):
    for item in win.items[:]:
        if not isinstance(item, Circle):
            item.undraw()


def main():
    f = open('results.json')
    data = json.load(f)

    win = GraphWin("Parallel Genetic TSP", data["width"], data["height"], autoflush=False)
    win.setBackground('white')

    win.getMouse()

    for city in data["points"]:
        draw_city(city[0], city[1], win)
    
    for generation in range(data["n_generations"]):
        clean_window(win)

        best_path = data["generations"][generation]["best_path"]
        for i in range(len(best_path)):
            city1 = data["points"] [best_path[i]]
            city2 = data["points"] [best_path[(i+1) % data["n_cities"]]]
            draw_line(city1[0], city1[1], city2[0], city2[1], win)
        
        draw_text(generation+1, data["generations"][generation]["best_path_distance"], win)

        update(60)
    
    win.getMouse()
    win.close()


if __name__ == "__main__":
    main()