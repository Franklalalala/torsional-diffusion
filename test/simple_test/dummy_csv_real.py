import csv


with open(r'test_smiles.csv') as rf, open(r'new_smiles.csv', 'a') as wf:
    csv_reader = csv.reader(rf)
    csv_writer = csv.writer(wf)
    for idx, a_line in enumerate(csv_reader):
        if idx > 300:
            break
        csv_writer.writerow(a_line)
        a = 1


# info = csv.reader(r'test_smiles.csv')
# for a_line in info:
#     a=1
a=1