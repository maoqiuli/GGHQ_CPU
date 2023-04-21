import os

def main():
    stage = "both"

    dataset = "spacev"
    datasize = 1
    attrsize = 16
    attrdim = 1
    M1 = 20
    M2 = 40 - M1


    command = "./main " + stage + " " + dataset + " " + str(datasize) + " " + str(attrsize) + " " + str(attrdim) + " " + str(M1) + " " + str(M2)

    os.system("cd build && make -j -B")
    os.system("cd build && " + command)


def all():
    os.system("cd build && make -j -B")

    stage = "both"
    datasetlist = ["spacev"]
    datasizelist = [1]
    attlist = [16]
    M1 = 20
    M2 = 40 - M1



    for i in range(len(datasetlist)):
        dataset = datasetlist[i]
        datasize = datasizelist[i]

        for attrsize in attlist:
            if attrsize == 12:
                attrdimlist = [1, 2, 4, 6, 8, 10]
            elif attrsize == 1024:
                attrdimlist = [1, 4, 16, 64, 128, 512]
            elif attrsize == 64:
                attrdimlist = [1, 2, 4, 8, 16, 32]
            elif attrsize == 16:
                attrdimlist = [1, 2, 4, 8, 16]
            else :
                attrdimlist = [1, 2, 4, 16, 64, 128]


            for attrdim in attrdimlist:
                if dataset == "sift":
                    process = "./main_char "
                else :
                    process = "./main "
                command = process + stage + " " + dataset + " " + str(datasize) + " " + str(attrsize) + " " + str(attrdim) + " " + str(M1) + " " + str(M2)
                os.system("cd build && " + command)


main()
# all()