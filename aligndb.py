import psycopg2
import numpy as np
db = psycopg2.connect(database='align170428')

def exe(sql, args=None):
    '''EXE - Execute a single SQL statement in its own transaction'''
    with db:
        with db.cursor() as c:
            c.execute(sql, args)

def nofail(sql, args=None):
    '''NOFAIL - Like EXE, but catches exceptions'''
    with db:
        with db.cursor() as c:
            try:
                c.execute(sql, args)
            except Exception as e:
                print(e)
                
def sel(sql, args=None):
    '''SEL - Run a SQL statement in its own transaction and return all
             fetched rows.'''
    with db:
        with db.cursor() as c:
            c.execute(sql, args)
            return c.fetchall()

X,Y,root = sel('select X,Y,root from info')[0]

def globaltotile(r, m, s, stage):
    '''GLOBALTOTILE - Return transformation matrix for given tile.
    t = GLOBALTOTILE(r, m, s, stage) returns a 3x3 affine matrix
    transforming global coordinates to pixel coordinates in the tile
    specified by (R, M, S) at the given processing stage. Processing
    stages are:
      BETA - Use only the translation from the BETAPOS table.
      UCT - Use the additional shift from the UCTSHIFT table.
    More stages will be defined as I progress.
    If you need the inverse transform, simply use np.linalg.inv.'''
    stage = stage.lower()
    t = np.eye(3)
    if stage=='beta':
        x, y = sel('''select xc, yc from betapos
                   where r=%s and m=%s and s=%s''', (r,m,s))[0]
        t[0,2] = X/2 - x
        t[1,2] = Y/2 - y
        return t
    elif stage=='uct':
        x, y = sel('''select xc, yc from betapos
                   where r=%s and m=%s and s=%s''', (r,m,s))[0]
        dx, dy = sel('''select dx, dy from uctshift
                   where r=%s and m=%s and s=%s''', (r,m,s))[0]
        t[0,2] = X/2 - x - dx
        t[1,2] = Y/2 - y - dy
        return t
    else:
        raise ValueError(f'Unknown stage: {stage}')

        
