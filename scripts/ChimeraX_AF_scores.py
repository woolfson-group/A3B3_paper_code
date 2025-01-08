# this script colours all atoms by plDDT scores in ChimeraX for all files
# loaded in a session used in Figure S12 (a,d,g)

from chimerax.core.commands import run

run(session, 'color bfactor palette 0,orange:49.999,orange:50,yellow:69.999,yellow:70,cornflowerblue:89.999,cornflowerblue:90,blue',)

for model in session.models:

    file_path = model.filename
    with open(file_path, 'r') as f:
        atom_lines = [line for line in f.readlines() if line.startswith('ATOM')]

    for line in atom_lines:
        line_split = line.split()
        atom_name = line_split[3]
        if atom_name == 'CA':

            chain = line_split[6]
            resnum = line_split[8]
            bfactor = float(line_split[-4])
            if bfactor <= 50:
                color = 'orange'
            elif 50 < bfactor <= 70:
                color = 'yellow'
            elif 70 < bfactor <= 90:
                color = 'cornflowerblue'
            else:
                color = 'blue'

            command = f'color #{model.id_string}/{chain}:{resnum}@{atom_name} {color} target c'
            run(session, command)
            