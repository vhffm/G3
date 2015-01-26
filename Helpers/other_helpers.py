def mkline(x1, y1, x2, y2):
    m = (y2 - y1) / (x2 - x1)
    n = y1 - m * x1
    return m, n
