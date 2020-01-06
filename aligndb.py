import psycopg2
import numpy as np

class DB:
    def __init__(self, **kwargs):
        if 'host' not in kwargs:
            kwargs['host'] = config.dbhost
        if 'user' not in kwargs:
            kwargs['user'] = config.dbuser
        if 'database' not in kwargs:
            kwargs['database'] = config.database
        self.db = psycopg2.connect(**kwargs)

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
        with self.db:
            with self.db.cursor() as c:
                c.execute(sql, args)
                ncols = len(c.description)
                raw = c.fetchall()
        raw = self.sel(sql, args)
        nrows = len(raw)
        if nrows==0:
            return [np.array([])] * ncols
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

    class Runlet:
        def __init__(self, r, s0, s1):
            self.r = r
            self.s0 = s0
            self.s1 = s1
        def __repr__(self):
            return f'Runlet(R{self.r} [{self.s0},{self.s1})'

    class RI:
        def nruns(self):
            '''NRUNS - Number of runs in the acquisition series
            NRUNS() returns the number of runs in the acquisition series.'''
            return self.R

        def nslices(self, r):
            '''NSLICES - Number of slices in a run
            NSLICES(r) returns the number of slices in given run.
            Note that runs count from one but slices within a run count from
            zero.'''
            return self.SS[r]

        def nmontages(self, r):
            '''NMONTAGES - Number of montages in a run
            NMONTAGES(r) returns the number of montages in given run.
            Note that runs count from one but montages within a run count from
            zero.'''
            return self.MM[r]

        def z0(self, r):
            '''Z0 - Z-coordinate of start of run
            Z0(r) returns the z-coordinate of the first slice of the given run.
            By definition, run one starts at z0 = 0.'''
            return self.zz0[r]

        def z(self, r, s):
            '''Z - Z-coordinate of a slice
            Z0(r, s) returns the z-coordinates of the given slice of the given
            run.'''
            z0 = self.zz0[r]
            return z0 + s

        def ncolumns(self, r):
            '''NCOLUMNS - Number of columns in a run
            NCOLUMNS(r) returns the number of columns in the given run.
            We use the config.py file to get this information.'''
            return config.ncolums(r)

        def nrows(self, r):
            '''NROWS - Number of rows in a run
            NROWS(r) returns the number of rows in the given run.
            Automatically determined from the number of columns.'''
            return self.nmontages(r) // self.ncolumns(r)

        def m2c(self, r, m):
            '''M2C - Get column number for a montage
            M2C(r, m) returns the column number for montage M in run R.'''
            return m % self.ncolumns(r)

        def m2r(self, r, m):
            '''M2C - Get row number for a montage
            M2C(r, m) returns the row number for montage M in run R.'''
            return m // self.ncolumns(r)

        def mright(self, r, m):
            '''MRIGHT - Find montage to the right of another
            MRIGHT(r, m) finds the montage number to the right of montage M
            in run R, or None if there is none.'''
            if self.m2c(r, m) < self.ncolumns(r)-1:
                return m+1
            else:
                return None

        def mleft(self, r, m):
            '''MLEFT - Find montage to the left of another
            MLEFT(r, m) finds the montage number to the left of montage M
            in run R, or None if there is none.'''
            if self.m2c(r, m) > 0:
                return m-1
            else:
                return None

        def mabove(self, r, m):
            '''MABOVE - Find montage above another
            MABOVE(r, m) finds the montage number above montage M
            in run R, or None if there is none. (Above means: in the
            negative y-direction.)'''
            if self.m2r(r, m) > 0:
                return m-self.ncolumns(r)
            else:
                return None

        def mbelow(self, r, m):
            '''MBELOW - Find montage below another
            MBELOW(r, m) finds the montage number below montage M
            in run R, or None if there is none. (Below means: in the
            positive y-direction.)'''
            if self.m2r(r, m) < self.nrows(r)-1:
                return m+self.ncolumns(r)
            else:
                return None
            
        def findz(self, z):
            '''Return a (r,s) pair for a given z'''
            for r0 in range(self.R):
                r = r0+1
                z0 = self.zz0[r]
                if z>=z0 and z<z0+self.SS[r]:
                    return (r, z-z0)
            return (None, None)

        def subvolume(self, z0, nz):
            '''Returns a list of runlets, spanning the given Z range'''
            rl = []
            r0,s0 = self.findz(z0)
            r1,s1 = self.findz(z0+nz-1)
            if r1 is None:
                return []
            s1 += 1
            for r in range(r0, r1+1):
                if r==r0:
                    s0a = s0
                else:
                    s0a = 0
                if r==r1:
                    s1a = s1
                else:
                    s1a = self.nslices(r)
                rl.append(DB.Runlet(r, s0a, s1a))
            return rl

    def runinfo(self):
        ri = self.sel('select r,M,S,z0 from runs')
        res = DB.RI()
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

