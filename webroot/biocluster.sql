/*
Navicat MySQL Data Transfer

Source Server         : 192.168.10.126
Source Server Version : 50173
Source Host           : localhost:3306
Source Database       : biocluster

Target Server Type    : MYSQL
Target Server Version : 50173
File Encoding         : 65001

Date: 2016-01-28 13:25:31
*/

SET FOREIGN_KEY_CHECKS=0;

-- ----------------------------
-- Table structure for `apilog`
-- ----------------------------
DROP TABLE IF EXISTS `apilog`;
CREATE TABLE `apilog` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `task_id` varchar(255) NOT NULL DEFAULT '0',
  `api` varchar(255) NOT NULL,
  `data` text,
  `addtime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `success` tinyint(4) NOT NULL DEFAULT '0',
  `last_update` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `server` varchar(255) DEFAULT NULL,
  `response` text,
  `response_code` varchar(255) DEFAULT NULL,
  `failed_times` int(11) NOT NULL DEFAULT '0',
  `reject` tinyint(4) NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  KEY `task_id_index` (`task_id`)
) ENGINE=MyISAM AUTO_INCREMENT=970 DEFAULT CHARSET=utf8;

-- ----------------------------
-- Records of apilog
-- ----------------------------

-- ----------------------------
-- Table structure for `clientkey`
-- ----------------------------
DROP TABLE IF EXISTS `clientkey`;
CREATE TABLE `clientkey` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `client` varchar(255) NOT NULL,
  `key` varchar(255) DEFAULT NULL,
  `ipmask` varchar(255) DEFAULT NULL,
  `timelimit` int(11) DEFAULT NULL,
  `max_workflow` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=4 DEFAULT CHARSET=utf8;

-- ----------------------------
-- Records of clientkey
-- ----------------------------
INSERT INTO `clientkey` VALUES ('1', 'client01', '1ZYw71APsQ', '172.16.3.0/24;192.168.10.0/24;127.0.0.1', '600', '500');
INSERT INTO `clientkey` VALUES ('2', 'test', 'Aw21cADS3u', '172.16.3.0/24;127.0.0.1', null, null);
INSERT INTO `clientkey` VALUES ('3', 'client02', '8A2q9C35Ts', '172.16.3.0/24;192.168.10.0/24', null, '5');

-- ----------------------------
-- Table structure for `pause`
-- ----------------------------
DROP TABLE IF EXISTS `pause`;
CREATE TABLE `pause` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `workflow_id` varchar(255) NOT NULL,
  `reason` varchar(255) DEFAULT '',
  `client` varchar(255) DEFAULT NULL,
  `ip` varchar(255) DEFAULT NULL,
  `add_time` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `has_pause` tinyint(4) NOT NULL DEFAULT '0',
  `pause_time` timestamp NULL DEFAULT NULL,
  `exit_pause` tinyint(4) NOT NULL DEFAULT '0',
  `exit_pause_time` timestamp NULL DEFAULT NULL,
  `has_continue` tinyint(4) NOT NULL DEFAULT '0',
  `continue_time` timestamp NULL DEFAULT NULL,
  `timeout` tinyint(4) NOT NULL DEFAULT '0',
  `timeout_time` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=MyISAM AUTO_INCREMENT=4 DEFAULT CHARSET=utf8;

-- ----------------------------
-- Records of pause
-- ----------------------------

-- ----------------------------
-- Table structure for `tostop`
-- ----------------------------
DROP TABLE IF EXISTS `tostop`;
CREATE TABLE `tostop` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `workflow_id` varchar(255) NOT NULL,
  `reson` varchar(255) NOT NULL DEFAULT '',
  `client` varchar(255) NOT NULL,
  `time` timestamp NULL DEFAULT CURRENT_TIMESTAMP,
  `ip` varchar(255) DEFAULT NULL,
  `done` tinyint(4) NOT NULL DEFAULT '0',
  `stoptime` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=11 DEFAULT CHARSET=utf8;

-- ----------------------------
-- Records of tostop
-- ----------------------------

-- ----------------------------
-- Table structure for `workflow`
-- ----------------------------
DROP TABLE IF EXISTS `workflow`;
CREATE TABLE `workflow` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `client` varchar(255) DEFAULT '' COMMENT '远程客户端名称',
  `workflow_id` varchar(255) NOT NULL COMMENT '流程ID',
  `json` text NOT NULL COMMENT '调用JSON文档',
  `ip` varchar(255) DEFAULT NULL COMMENT '远程客户端IP地址',
  `add_time` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `has_run` tinyint(4) NOT NULL DEFAULT '0',
  `run_time` timestamp NULL DEFAULT NULL,
  `is_end` tinyint(4) NOT NULL DEFAULT '0',
  `end_time` timestamp NULL DEFAULT NULL,
  `server` varchar(255) DEFAULT NULL,
  `last_update` timestamp NULL DEFAULT NULL,
  `is_error` tinyint(4) DEFAULT '0',
  `error` varchar(255) DEFAULT NULL,
  `output` varchar(255) DEFAULT NULL,
  `pid` int(11) DEFAULT NULL,
  `workdir` varchar(255) DEFAULT NULL,
  `paused` tinyint(4) NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  UNIQUE KEY `workflow_id_index` (`workflow_id`) USING BTREE
) ENGINE=InnoDB AUTO_INCREMENT=336 DEFAULT CHARSET=utf8;

-- ----------------------------
-- Records of workflow
-- ----------------------------
