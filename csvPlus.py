import csv

class DictReader():
    """csvPlus.DictReader is a wrapper for csv.DictReader, with the additional method reset"""

    def __init__(self, datafile):
        self.rawdata = datafile.readlines()
        self.data = csv.DictReader(self.rawdata)
        self.index = 0

    def __iter__(self):
        """wrap the __iter__() method of csv.DictReader"""
        return self.data

    def next(self):
        """wrap the next() method of csv.DictReader"""
        return self.data.next()

    def reset(self):
        """reset the DictReader iterator to the beginning of the data file"""
        self.data = None
        self.data = csv.DictReader(self.rawdata)
        return True

class reader():
    """csvPlus.reader is a wrapper for csv.DictReader, with the additional method reset"""

    def __init__(self, datafile):
        self.rawdata = datafile.readlines()
        self.data = csv.reader(self.rawdata)
        self.index = 0

    def __iter__(self):
        """wrap the __iter__() method of csv.DictReader"""
        return self.data

    def next(self):
        """wrap the next() method of csv.DictReader"""
        return self.data.next()

    def reset(self):
        """reset the DictReader iterator to the beginning of the data file"""
        self.data = None
        self.data = csv.reader(self.rawdata)
        return True
