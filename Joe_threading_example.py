import threading
import Queue
import log
import traceback


class Worker:
    def __init__(self, numthreads):
        self.numthreads = numthreads
        self.running = False
        self.q = Queue.Queue()
        self.threads = []

    def clear(self):
        try:
            while not self.q.empty():
                self.q.get(False)
        except Exception, e:
            pass
        
    def start(self):
        self.running = True
        for i in xrange(self.numthreads):
            t = threading.Thread(target=self.runloop)
            t.setDaemon(True)
            t.start()
            self.threads.append(t)

    def stop(self):
        self.running = False
        q.join()
        for t in self.threads:
            t.join()

    def runloop(self):
        while (self.running):
            item = self.q.get()
            try:
                item()
            except:
                log.warn(traceback.format_exc())
            self.q.task_done()

    def add(self, func):
        self.q.put(func)
