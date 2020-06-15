from ROOT import *

canvas = TCanvas("Canvas_1", "Canvas_1", 533, 76, 1383, 852)


h1 = TH1F("edep", "edep", 390, 0, 0.39)

with open("out/photons.txt", "r") as f:
    lines = f.readlines()
    for line in lines:
        data = line.split(" ")
        h1.Fill(float(data[3]))


h1.Draw()
canvas.SaveAs("edep.png")
