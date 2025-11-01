from FileRead import readcol
import sys

[star_id, x_sw, y_sw, x_lw, y_lw] = readcol(sys.argv[1], 'f,ffff')

for sw_lw in ["sw", "lw"]:
    f = open("ds9_" + sys.argv[1].split(".")[0] + "_" + sw_lw + ".reg", 'w')
    f.write("""# Region file format: DS9 version 4.1
    global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
    image
    """)

    for i in range(len(star_id)):
        if sw_lw == "sw":
            f.write("circle(%i,%i,10)\n" % (x_sw[i], y_sw[i]))
        else:
            f.write("circle(%i,%i,10)\n" % (x_lw[i], y_lw[i]))

    f.close()
