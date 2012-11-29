class int_incrementer:
    def __init__(self, startval):
        self.current = hash(startval)

    def inc(self):
        self.current += 1
        return self.current
    
counter = int_incrementer(0)
