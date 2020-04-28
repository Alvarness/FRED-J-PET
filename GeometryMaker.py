
# scin_dim_x = 2.5
# scin_dim_y = 0.6
# scin_dim_z = 50
# space_between = 0.1

# scin_in_module = 13
# modules = 24
# radius = 38.186

# f = open("JPetModular.txt", "w")

# f.write("material< \n\
#     ID = ej230 \n\
#     rho = 1.023\n\
#     Ipot = 9999\n\
#     composition = [C,H]\n\
#     fractions = [47.619,52.381]\n\
# material>\n \n \n")

# for module in range(modules):
#     group = "group: module_{} ".format(module)
#     for scin in range(scin_in_module):
#         f.write("region: scintillator_{}_{} ; pivot=[0.5,0.5,0.5] ; L=[{},{},{}] ; O=[{}, {}, 0] ; material = water\n".format(
#             module, scin, scin_dim_x, scin_dim_y, scin_dim_z, radius, scin * 0.7 - 0.7 * 6))
#         group += "scintillator_{}_{} ".format(module, scin)
#     group += "\n"
#     f.write(group)
#     f.write("transform: module_{} rotate z {} \n".format(
#         module, module * 360.0 / modules))

# f.close()

scin_dim_x = 1.9
scin_dim_y = 0.7
scin_dim_z = 50

layer1 = 48
layer2 = 48
layer3 = 96
radius1 = 42.5
radius2 = 46.75
radius3 = 57.5


f = open("JPetBigBarrel.txt", "w")

for scintillator in range(layer1):
    f.write(
        "region: scintillator_1_{} ; pivot=[0.5,0.5,0.5] ; L=[{},{},{}] ; O=[{}, 0, 0] ; material = water\n".format(
            scintillator, scin_dim_x, scin_dim_y, scin_dim_z, radius1))
    f.write("transform: scintillator_1_{} rotate z {} \n".format(
        scintillator, scintillator * (360.0 / layer1)))

for scintillator in range(layer2):
    f.write(
        "region: scintillator_2_{} ; pivot=[0.5,0.5,0.5] ; L=[{},{},{}] ; O=[{}, 0, 0] ; material = water\n".format(
            scintillator, scin_dim_x, scin_dim_y, scin_dim_z, radius2))
    f.write("transform: scintillator_2_{} rotate z {} \n".format(
        scintillator, scintillator * (360.0 / layer2) + 3.75))


for scintillator in range(layer3):
    f.write(
        "region: scintillator_3_{} ; pivot=[0.5,0.5,0.5] ; L=[{},{},{}] ; O=[{}, 0, 0] ; material = water\n".format(
            scintillator, scin_dim_x, scin_dim_y, scin_dim_z, radius3))
    f.write("transform: scintillator_3_{} rotate z {} \n".format(
        scintillator, scintillator * (360.0 / layer3) + 1.875))


f.close()
