<?php
function checkIdentity($code, $fileName, $mode){
	$info = array();
	if ($mode == "tsanger"){
	    $servername = "192.168.10.51";
	    $username = "i-sanger";
	    $password = "sg123123";
	}
	else if ($mode == "sanger"){
	    $servername = "192.168.10.179";
	    $username = "tpuser";
		$password = "uDB#345TT";
	}
	else {
	    $info["success"] = false;
		$info["info"] = "模式错误";
		return $info;
	}
	if ($mode == "tsanger"){
	    $database = "tp";
	}
	else{
	    $database = "isanger";
	}
	$conn = new mysqli($servername, $username, $password, $database);
	if ($conn->connect_error) {
	    die("连接失败: " . $conn->connect_error);
	}
	echo $conn->connect_error;
	$stmt = $conn->prepare("SELECT create_time, related_task_id FROM sg_download_code WHERE code=?");
	$stmt->bind_param("s", $code);
	$stmt->execute();
	$stmt->bind_result($create_time, $task_id);
	$result = array();
	while($stmt->fetch()){
		$result[$task_id] = $create_time;
		$tmp_id = $task_id;
		$tmp_time = $create_time;
		$tmp_time = DateTime::createFromFormat('Y-m-d H:i:s', $tmp_time);
	}
	if (count($result) == 0){
		$info["success"] = false;
		$info["info"] = "验证码" . $code . "错误";
	    return $info;
	}
	else{
	    $time_now = new DateTime();
		if ($tmp_time->add(new DateInterval("PT15H"))< $time_now){
		    $info["success"] = false;
			$info["info"] = "验证码已超时";
			return $info;
		}
		else{
			if (!preg_match("/" . $tmp_id . "/", $fileName)){
				$info["success"] = false;
				$info["info"] = "文件路径不正确";
				return $info;
			}
			else{
			    $info["success"] = true;
				$info["info"] = "验证通过";
				return $info;
			}
		}
	}
}



$file = $_POST['file'];
$code = $_POST['indentity_code'];
$mode = $_POST['mode'];

$filename = basename($file);
$encoded_filename = urlencode($filename);

$ua = $_SERVER["HTTP_USER_AGENT"];
if ((preg_match("/MSIE/", $ua)) or (preg_match("/Firefox/", $ua)) or (preg_match("/Chrome/", $ua))){
	throw new Exception("不可以使用浏览器访问！");
}
else {
    $info = checkIdentity($code, $file, $mode);
	if (!$info["success"]){
		header('HTTP/1.1 405 '. $info["info"]);
	}
	else{
        header('Content-Type: application/octet-stream');
	    header('Content-Disposition: attachment; filename="' . $filename . '"');
	}
}

header("X-Sendfile: $file");
exit();
?>
