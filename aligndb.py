import psycopg2
import numpy as np

class DB:
    def __init__(self, **kwargs):
        if 'host' not in kwargs:
            kwargs['host'] = 'localhost'
        self.db = psycopg2.connect(database='align170428', **kwargs)
        self.X, self.Y, self.root = self.sel('select X,Y,root from info')[0]

    def exe(self, sql, args=None):
        '''EXE - Execute a single SQL statement in its own transaction'''
        with self.db:
            with self.db.cursor() as c:
                c.execute(sql, args)
    
    def nofail(self, sql, args=None):
        '''NOFAIL - Like EXE, but catches exceptions'''
        with self.db:
            with self.db.cursor() as c:
                try:
                    c.execute(sql, args)
                except Exception as e:
                    print(e)
                    
    def sel(self, sql, args=None):
        '''SEL - Run a SQL statement in its own transaction and return all
                 fetched rows.'''
        with self.db:
            with self.db.cursor() as c:
                c.execute(sql, args)
                return c.fetchall()
    
    def vsel(self, sql, args=None):
        '''VSEL - Like SEL, but returns results as a list of numpy arrays.
        VSEL(sql, args) executes the SQL code, interpolating ARGS for "%s"
        placeholders, and returns the results as a list of numpy arrays.'''
        raw = self.sel(sql, args)
        nrows = len(raw)
        if nrows==0:
            return None
        ncols = len(raw[0])
        res = []
        for c in range(ncols):
            if type(raw[0][c])==int:
                ar = np.zeros(nrows, dtype=int)
            else:
                ar = np.zeros(nrows)
            for r in range(nrows):
                ar[r] = raw[r][c]
            res.append(ar)
        return res
   
    def globaltotile(self, r, m, s, stage):
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
            x, y = self.sel('''select xc, yc from betapos
                       where r=%s and m=%s and s=%s''', (r,m,s))[0]
            t[0,2] = X/2 - x
            t[1,2] = Y/2 - y
            return t
        elif stage=='uct':
            x, y = self.sel('''select xc, yc from betapos
                       where r=%s and m=%s and s=%s''', (r,m,s))[0]
            dx, dy = self.sel('''select dx, dy from uctshift
                       where r=%s and m=%s and s=%s''', (r,m,s))[0]
            t[0,2] = X/2 - x - dx
            t[1,2] = Y/2 - y - dy
            return t
        else:
            raise ValueError(f'Unknown stage: {stage}')

    def runinfo(self):
        ri = self.sel('select r,M,S,z0 from runs')
        class RI:
            def nruns(self):
                return self.R
            def nslices(self, r):
                return self.SS[r]
            def nmontages(self, r):
                return self.MM[r]
            def z0(self, r):
                return self.zz0[r]
            def ncolumns(self, r):
                if self.MM[r]<=3:
                    return 1
                else:
                    return 2
            def nrows(self, r):
                return self.nmontages(r) // self.ncolumns(r)
        res = RI()
        res.R = len(ri)
        res.MM = {}
        res.SS = {}
        res.zz0 = {}
        for row in ri:
            r = row[0]
            res.MM[r] = row[1]
            res.SS[r] = row[2]
            res.zz0[r] = row[3]
        return res
    
