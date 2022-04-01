# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import pandas as pd
import json
import time

class TfNetworkWorkflow(Workflow):
    """
    TfNetwork description
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TfNetworkWorkflow, self).__init__(wsheet_object)
        options = list()
        for each in self._sheet.options():
            options.append(dict(name=each, type="string"))
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.tf_network = self.add_tool("ref_rna.tf_network")
        self.dump_tool = self.api.api("lnc_rna.api_base")

    def run(self):
        self.start_listener(); self.fire("start") # if you have no tools, you should use this line
        self.run_tf_network()
        self.set_db()
        # super(TfNetworkWorkflow, self).run()

    def run_tf_network(self):
        time.sleep(7)
        tfbs_predict = pd.read_table(self.option("tfbs_predict"), header=0)
        tf_network = tfbs_predict.drop_duplicates(["tf_geneid", "gene_name_tf", "target_id", "gene_name_target"])
        self2other_ind = tf_network['tf_geneid'] != tf_network["target_id"]
        tf_network = tf_network[self2other_ind]
        if "corr" in tf_network.columns:
            target_cols = ["tf_geneid", "gene_name_tf", "target_id", "gene_name_target",
                           "corr", "corr_pvalue", "corr_padjust", "p-value", "q-value"]
        else:
            target_cols = ["tf_geneid", "gene_name_tf", "target_id", "gene_name_target", "p-value", "q-value"]
        tf_network = tf_network.loc[:, target_cols]
        # filtering
        if self.option("qv_thresh") == "0":
            tf_network = tf_network[tf_network['p-value'] <= float(self.option("thresh"))]
        else:
            tf_network = tf_network[tf_network['q-value'] <= float(self.option("thresh"))]
        if 'corr' not in tf_network.columns:
            tf_network = tf_network[tf_network['corr'].abs() >= float(self.option("corr_cutoff"))]
            if self.option("corr_use_pdajust") == "0":
                tf_network = tf_network[tf_network['corr_pvalue'] <= self.option("corr_pvalue")]
            else:
                tf_network = tf_network[tf_network['corr_padjust'] <= self.option("corr_pvalue")]
        tf_network.to_csv(self.output_dir + '/network.edges.txt', sep='\t', index=False, header=True)
        tf_nodes = tf_network.loc[:, ["tf_geneid", "gene_name_tf"]].drop_duplicates()
        tf_nodes.columns = ['id', 'name']
        tf_nodes['group'] = 'TF'
        gene_nodes = tf_network.loc[:, ["target_id", "gene_name_target"]].drop_duplicates()
        gene_nodes.columns = ['id', 'name']
        gene_nodes['group'] = 'Gene'
        all_nodes = tf_nodes.append(gene_nodes)
        duplicated_ind = all_nodes.duplicated(['id', 'name'])
        all_nodes['group'][duplicated_ind] = "Gene_TF"
        dup_nodes = all_nodes[duplicated_ind].drop_duplicates()
        all_nodes = all_nodes.drop_duplicates(['id', 'name'], keep=False)
        all_nodes = all_nodes.append(dup_nodes)
        all_nodes.to_csv(self.output_dir + '/network.nodes.txt', header=True, index=False, sep='\t')
        with open(self.output_dir + "/network.json", "w") as f:
            tmp_nodes = all_nodes.set_index("id")
            tmp_nodes['index'] = list(range(0, tmp_nodes.shape[0]))
            index_node_pd = tmp_nodes.loc[:, ["index"]]
            tmp_links = tf_network.loc[:, ['tf_geneid', 'target_id']]
            if 'corr' in tf_network.columns:
                tmp_links['weight'] = tf_network['corr']
            else:
                tmp_links['weight'] = 1
            tmp_links.columns = ["source", "target", "distance"]
            tmp_links['distance'] = tmp_links['distance'].abs()
            tmp_links = tmp_links.join(index_node_pd, on="source")
            tmp_links = tmp_links.join(index_node_pd, on="target", lsuffix="_source", rsuffix="_target")
            tmp_links = tmp_links.loc[:, ["index_source", "index_target", "distance"]]
            tmp_links.columns = ["source", "target", "distance"]
            json.dump(dict(
                links=json.loads(tmp_links.to_json(orient="records")),
                nodes=all_nodes.to_dict("records")), f, indent=2,
            )

    def set_db(self):
        """
        dump data to db
        """
        workflow_output = self.get_workflow_output_dir()
        # add result info
        # self.dump_tool.add_tf_network_detail(self.tf_network.output_dir, self.option("main_id"))
        self.dump_tool.update_db_record('sg_tf_network', self.option('main_id'), output_dir=workflow_output, status="end")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "TfNetwork", 0, "211242"],
        ])
        super(TfNetworkWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output
