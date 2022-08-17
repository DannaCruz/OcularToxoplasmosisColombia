import csv
import pathlib


def processFile(name):
    f = open('./BaseFiles/' + name + '.data', "r")
    header = f.readline()
    line = f.readline().split('\n')[0]
    rows = []
    while (line != ''):
        partes = line.split('|')
        fecha = partes[0].split(' ')[0]
        hora = partes[0].split(' ')[1]
        precipitacion = partes[1].split('\n')[0]
        row = [fecha, hora, float(precipitacion)]
        if'2015'in fecha or '2016'in fecha or '2017'in fecha or'2018'in fecha or'2019'in fecha:
            rows.append(row)
        line = f.readline().split('\n')[0]
    nheader = [header.split('|')[0], 'Hora', header.split('|')[1].split('\n')[0]]
    if len(rows)>0:
        d = open('./ProcessFiles/' + name + '.csv', 'w')
        writer = csv.writer(d)
        writer.writerow(nheader)
        writer.writerows(rows)
        d.close()
    f.close()

def openAllFiles(directory):
    i = 1
    for path in pathlib.Path(directory).iterdir():
        if path.is_file():
            processFile(path.__str__().split('.')[0].split('\\')[1])
            print('File number', i, 'processed')
            i += 1

if __name__ == '__main__':
    openAllFiles('./BaseFiles')

