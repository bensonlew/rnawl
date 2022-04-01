# -*- coding: utf-8 -*-

from __future__ import print_function

if __name__ == '__main__':
    states_needed = {'mt', 'wa', 'or', 'id', 'nv', 'ut', 'ca', 'az'}
    states_dict = {'kone': {'id', 'nv', 'ut'},
                   'ktwo': {'wa', 'id', 'mt'},
                   'kthree': {'or', 'nv', 'ca'},
                   'kfour': {'nv', 'ut'},
                   'kfive': {'ca', 'az'}}
    final_stations = set()
    while states_needed:
        best_station = None
        states_covered = set()
        for station, states_for_station in states_dict.items():
            covered = states_needed & states_for_station
            if len(covered) > len(states_covered):
                best_station = station
                states_covered = covered
        final_stations.add(best_station)
        states_needed -= states_covered
    print(final_stations)
