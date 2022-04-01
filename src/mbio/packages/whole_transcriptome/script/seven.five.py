# -*- coding: utf-8 -*-

from __future__ import print_function


def find_lowest_cost_node(cost_dict, processed):
    _dict = {k: v for k, v in cost_dict.items() if k not in processed}
    if _dict:
        return min(_dict, key=lambda k: _dict[k])


def complete_parent(graph, cost_dict, parent_dict):
    processed = {'end'}
    node = find_lowest_cost_node(cost_dict, processed)
    while node:
        cost = cost_dict[node]
        neighbor_dict = graph[node]
        for n, c in neighbor_dict.items():
            new_cost = cost + c
            if cost_dict[n] > new_cost:
                cost_dict[n] = new_cost
                parent_dict[n] = node
        processed.add(node)
        node = find_lowest_cost_node(cost_dict, processed)


def display_route(parent_dict):
    route = ['end']
    node = parent_dict['end']
    while node != 'start':
        route.append(node)
        node = parent_dict[node]
    route.append('start')
    for i, n in enumerate(route[::-1]):
        if i != len(route) - 1:
            print('{} ->'.format(n), end=' ')
        else:
            print(n)


if __name__ == '__main__':
    graph = {'start': {'A': 6, 'B': 2}, 'A': {'end': 1}, 'B': {'A': 3, 'end': 5}}
    cost_dict = {'A': 6, 'B': 2, 'end': float('inf')}
    parent_dict = {'A': 'start', 'B': 'start', 'end': None}
    complete_parent(graph, cost_dict, parent_dict)
    display_route(parent_dict)
