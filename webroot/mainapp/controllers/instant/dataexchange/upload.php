<html>
<head>
<title>文件上传</title>
</head>
<body>
<?php
// function checkIdentity($code, $input_dir, $mode){
// 	$info = array();
// 	if ($mode == "tsanger"){
// 	    $servername = "192.168.10.51";
// 	    $username = "i-sanger";
// 	    $password = "sg123123";
// 	}
// 	else if ($mode == "sanger"){
// 		$servername = "192.168.10.179";
// 		$username = "tpuser";
// 		$password = "uDB#345TT";
// 	}
// 	else{
// 	    $info["success"] = false;
// 		$info["info"] = "模式错误";
// 		return $info;
// 	}
//     $database = "tp";
//     $conn = new mysqli($servername, $username, $password, $database);
//     if ($conn->connect_error) {
// 		die("连接失败: " . $conn->connect_error);
//     }
//     echo $conn->connect_error;
// 	$stmt = $conn->prepare("SELECT create_time, rel_dir FROM sg_upload_code WHERE code=?");
// 	$stmt->bind_param("s", $code);
// 	$stmt->execute();
// 	$stmt->bind_result($create_time, $rel_dir);
//     $result = array();
// 	while($stmt->fetch()){
// 		array_push($result, $rel_dir);
// 	    $tmp_time = $create_time;
// 		$tmp_time = DateTime::createFromFormat('Y-m-d H:i:s', $tmp_time);
// 	}
// 	if (count($result) == 0){
// 		$info["success"] = false;
// 		$info["info"] = "验证码" . $code . "错误";
// 		return $info;
// 	}
// 	else{
// 	    $time_now = new DateTime();
// 		if ($tmp_time->add(new DateInterval("PT15H"))< $time_now){
// 			$info["success"] = false;
// 			$info["info"] = "验证码已超时";
// 			return $info;
// 		}
// 		else{
// 			if(!in_array($input_dir, $result)){
// 				$info["success"] = false;
// 				$info["info"] = "文件路径不正确";
// 				return $info;
// 			}
// 			else{
// 		        $info["success"] = true;
// 			    $info["info"] = "验证通过";
// 			    return $info;
// 			}
// 		}
// 	}
// }

var_dump($_FILES);
var_dump($_POST);
// $code = $_POST["code"];
// $mode = $_POST["mode"];
$target_path = $_POST["target_path"];
$input_dir = $_POST["rel_dir"];
$dir_name = dirname($target_path);

$verify = $_POST['verify'];
if(empty($verify) || md5($verify) != 'ce85e3a78883451d933aeaa2b7f38d76') {
	exit('非法访问');
}

$ua = $_SERVER["HTTP_USER_AGENT"];
if ((preg_match("/MSIE/", $ua)) or (preg_match("/Firefox/", $ua)) or (preg_match("/Chrome/", $ua))){
    throw new Exception("不可以使用浏览器访问！");
}
else{
	// $info = checkIdentity($code, $input_dir, $mode);
	// if (!$info["success"]){
	//     header('HTTP/1.1 405 '. $info);
	// }
	// else{
		if (!is_dir($dir_name)){
		    mkdir($dir_name, 0777, true);
		}
		move_uploaded_file($_FILES["sources"]["tmp_name"], $target_path);
		chmod($target_path, 0777);
		echo "上传完成！";
	// }
}
?>
</body>
</html>
