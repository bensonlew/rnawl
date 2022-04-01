# -*- coding: utf-8 -*-
# __author__ = 'kefei.huang'

dbconfigure={'host':'10.100.200.135','name':'mjlimsbak','passwd':'Q6k9pU9ZN5FTow','database':'mjlims'}
DicWQ = {}
DicWS = {}
WQ_sample_name = []
WQ_error_count = 0
with codecs.open("/mnt/ilustre/users/sanger-dev/workspace/20171109/Single_file_check_test/FileCheckV4/20171019.split.csv","r","utf-8") as IN:
	for line in IN:
		linetemp = line.split(",")
		if linetemp[0].startswith("WQ"):
			if re.search(r'-S|-M|-F',linetemp[0]) == None: ## -s -f -m
				WQ_error_count += 1
			case_name = linetemp[0].split("-")[0]
			DicWQ[case_name] = 1
			WQ_sample_name.append(linetemp[0])
		if linetemp[0].startswith("WS"):
			DicWS[linetemp[0]] = 1

###链接数据库
con = mdb.connect(dbconfigure["host"], dbconfigure["name"], dbconfigure["passwd"], dbconfigure["database"], charset = 'utf8')
cur = con.cursor()

WQ_sample = "','".join(DicWQ.keys())
WS_sample = "','".join(DicWS.keys())
sql = "SELECT i.note, o.create_date, o.name, i.patient_name, st.name, i.code, i.parented, o.sjrq, " \
		"i.accept_date, o.isurgent, o.gestational_weeks, t.name, t.principal FROM sample_info i " \
		"LEFT JOIN sample_order o ON i.sample_order = o.id LEFT JOIN primary_task t ON t.id = o.advance " \
		"LEFT JOIN dic_sample_type st ON st.id = i.sample_kind WHERE o.id IS NOT NULL AND i.note LIKE 'WQ%'" \
		"and i.note in ('{}') ORDER BY o.create_date DESC".format(WQ_sample)
#WQ17093590	2017-09-20 14:25:37	张瑶瑶	黄小明	精斑	JB1709200003	F1	2017-09-16 14:30:23	2017-09-20 14:24:57	0	26	曾蓉（网络客服）	曾蓉  示例
cur.execute(sql)
case_infor = cur.fetchall()

client = pymongo.MongoClient('mongodb://192.168.10.189:27017/')
db = client['tsanger_paternity_test_v2']
collection = db['sg_update_test']

for each in case_infor:
	insert_data = {
	  "case_name": each[0],
	  "create_time": each[1],
	  "ask_person": each[2],
	  "sample_name": each[3],
	  "sample_type": each[4],
	  "sample_number": each[5],
	  "sample_id": each[6],
	  "ask_time": each[7],
	  "accept_time": each[8],
	  "gestation_week": each[10],
	  "company_name": each[11],
	  "contacts_people": each[12],
	  "name" : each[0] + "-" + each[6],
	  'insert_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	}
	status = scollection.update({"name":each[0] + "-" + each[6]},{'$set':insert_data},upsert=True)