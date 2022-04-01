# -*- coding: utf-8 -*-

from collections import deque

person_is_seller = lambda name: name[-1] == 'm'


def has_mango_seller(search_queue, graph):
    searched = set()
    while search_queue:
        person = search_queue.popleft()
        if person in searched:
            continue
        else:
            searched.add(person)
        if person_is_seller(person):
            print('{} is a mango seller!'.format(person))
            return True
        else:
            search_queue.extend(graph[person])
    return False


if __name__ == '__main__':
    graph = dict()
    graph["you"] = ["alice", "bob", "claire"]
    graph["bob"] = ["anuj", "peggy"]
    graph["alice"] = ["peggy"]
    graph["claire"] = ["thom", "jonny"]
    graph["anuj"] = []
    graph["peggy"] = []
    graph["thom"] = []
    graph["jonny"] = []

    search_queue = deque()
    search_queue.extend(graph['you'])

    has_mango_seller(search_queue, graph)
