import os
import json
import shutil

from commonFunction import re_path

''' --------------------- ДОПОЛНИТЕЛЬНЫЕ ФУНКЦИИ ДЛЯ ЛИБЫ --------------------
ПЕРЕСМОТРЕТЬ И ПЕРЕПИСАТЬ

def recursive_copy(src, dst, overwrite=False):
    # modified https://www.reddit.com/r/learnpython/comments/6g67o8/how_to_recursively_copy_directory_contents_into/
    #    os.chdir(src)
    src = re_path(src)
    dst = re_path(dst)
    if not os.path.exists(dst):
        os.makedirs(dst)

    for item in os.listdir(src):
        full_src = os.path.join(src, item)
        full_dst = os.path.join(dst, item)

        if os.path.isfile(full_src):
            if overwrite and os.path.isfile(full_dst):
                os.remove(full_dst)
            shutil.copy2(full_src, full_dst)

        elif os.path.isdir(full_src):
            if not os.path.exists(full_dst):
                os.makedirs(full_dst)
            recursive_copy(full_src, full_dst, overwrite=overwrite)


def recursive_move(src, dst, overwrite=False):

    src = re_path(src)
    dst = re_path(dst)
    if not os.path.exists(dst):
        os.makedirs(dst)

    for item in os.listdir(src):
        full_src = os.path.join(src, item)
        full_dst = os.path.join(dst, item)

        if os.path.isfile(full_src):
            if overwrite and os.path.isfile(full_dst):
                os.remove(full_dst)

            shutil.move(full_src, full_dst)

        elif os.path.isdir(full_src):
            if not os.path.exists(full_dst):
                os.makedirs(full_dst)
            else:
                shutil.rmtree(full_dst)

            recursive_move(full_src, full_dst)
            shutil.rmtree(full_src)

def copy_(self, src='', dst='', rename='', putInto='0'):
    """
     Копировать файл или папку.

    :param src:      - путь, указывающий на файл или папку, которые надо
                       скопировать.
    :param dst:      - путь, указывающий куда копировать.
    :param rename:   - переименовать файл src. Если rename="new_name.xml",
                       то это новое имя объекта src.
    :param putInto: - функция, определяющая выполнение команды copy. При-
                       нимает значения True or False, '1' или '0'.
                       Подробности в спецификации на git.
    """

    if os.path.isdir(src):
        if os.path.isdir(dst) and putInto == '0':
            recursive_copy(src, dst)

        elif os.path.isdir(dst) and putInto == '1':
            dst = re_path(os.path.join(DST, os.path.split(SRC)[-1]))
            recursive_copy(src, dst)

    elif os.path.isfile(src):
        if os.path.isfile(dst) and putInto == '0':
            shutil.copy2(src, dst)

        elif os.path.isdir(dst) and putInto == '1' and rename:
            dst = re_path(os.path.join(dst, rename))
            shutil.copy2(src, dst)


def move_(self, src='', dst='', rename='', putInto='0'):
    """
     Переместить файл или папку.

    :param src:      - путь, указывающий на файл или папку, которые надо
                       переместить.
    :param dst:      - путь, указывающий куда переместить.
    :param rename:   - переименовать файл src. Если rename="new_name.xml",
                       то это новое имя объекта src.
    :param putInto: - функция, определяющая выполнение команды move. При-
                       нимает значения True or False, '1' или '0'.
                       Подробности в спецификации на git.
    """

    if os.path.isdir(src):
        if os.path.isdir(dst) and putInto == '0':
            recursive_move(src, dst)

        elif os.path.isdir(dst) and putInto == '1':
            dst = re_path(os.path.join(DST, os.path.split(SRC)[-1]))
            recursive_move(src, dst)
            shutil.rmtree(src)

    elif os.path.isfile(src):
        if os.path.isfile(dst) and putInto == '0':
            shutil.move(src, dst)

        elif os.path.isdir(dst) and putInto == '1' and rename:
            dst = re_path(os.path.join(dst, rename))
            shutil.move(src, dst)

    # else:
    #     self.taskMsg.append(
    #         'test "' + self.__class__.taskNamesList[
    #             self.__class__.taskIndex] + '" data moving error: "' +
    #         fullPath[0] + '" missed, no files moved')
    #     self.testErrCodes.append(1)


def delete_(self, path='', include='', _except=''):
    """
     Удалить файл(-ы) или папку(-и). Удалить всё с расширением (*.__),
     удалить все кроме "name_file.ext" или "*.__"

    :param path:     - путь к файлу или папке.
    :param include:  - Параметр, в котором пишется, что надо удалить.
                         Используется, когда path указывает на директорию.
    :param _except:  - имена/тип файлов, которые не удалять.
                                                    (_except='*.wav')
    """

    if os.path.exists(path):
        if include:
            if include.startswith('*.') or include.startswith('.'):

                list_of_files_to_del_by_extension = \
                    [item.strip() for item in include.split(',')]

                for ext in list_of_files_to_del_by_extension:
                    for file in os.listdir(path):
                        full_path = re_path(os.path.join(path, file))

                        if file.endswith(ext.strip('*')) and \
                                os.path.isfile(full_path):
                            os.remove(full_path)

            else:

                list_of_files_or_dir_to_del = \
                    [item.strip() for item in include.split(',')]

                for file in list_of_files_or_dir_to_del:
                    full_path = re_path(os.path.join(path, file))

                    if os.path.isfile(full_path):
                        os.remove(full_path)

                    elif os.path.isdir(full_path):
                        shutil.rmtree(full_path)

        elif _except:

            list_of_files_to_del = \
                [item.strip() for item in _except.split(',')]

            for file in os.listdir(path):
                for item in list_of_files_to_del:
                    full_path = re_path(os.path.join(path, item))

                    if (item.startswith('*.') or item.startswith('.')):
                        ext_file = os.path.splitext(file)[-1]
                        list_ext = [
                            ext.strip('*') for ext in list_of_files_to_del
                        ]
                        if ext_file in list_ext:
                            continue
                        else:
                            self.delete_(path=path, include=file)

                    elif (os.path.isfile(full_path) or
                          os.path.isdir(full_path)) and \
                            file not in list_of_files_to_del:
                        self.delete_(path=path, include=file)

        else:
            shutil.rmtree(path)

    else:
        print('Сообщение об ошибке')
        return

    # else:
    #     self.taskMsg.append(
    #         'test "' + self.__class__.taskNamesList[
    #             self.__class__.taskIndex] + '" data prepareing error: "'
    # + path + '" missed, no files deleted')
    #     self.testErrCodes.append(1)
    #     return


def clear_(self, path='', _except=''):
    """
     Очистить директорию, очистить все кроме файла, файлов.

    :param path:    - путь к папке.
    :param _except: - имена/тип файлов, которые не удалять.
                                                    (_except='*.wav')
    """

    if os.path.exists(path):
        if _except:

            list_of_files_to_del = \
                [item.strip() for item in _except.split(',')]

            for file in os.listdir(path):
                for item in list_of_files_to_del:
                    full_path = re_path(os.path.join(path, item))

                    if (item.startswith('*.') or item.startswith('.')):
                        ext_file = os.path.splitext(file)[-1]
                        list_ext = [
                            ext.strip('*') for ext in list_of_files_to_del
                        ]
                        if ext_file in list_ext:
                            continue
                        else:
                            self.delete_(path=path, include=file)

                    elif (os.path.isfile(full_path) or
                          os.path.isdir(full_path)) and \
                            file not in list_of_files_to_del:
                        self.delete_(path=path, include=file)

        else:
            for f_name in os.listdir(path):
                full_path = os.path.join(path, f_name)

                if os.path.isfile(full_path):
                    os.remove(full_path)
                elif os.path.isdir(full_path):
                    shutil.rmtree(full_path)

    else:
        print('Сообщение об ошибке')
        return

    # else:
    #     self.taskMsg.append(
    #         'test "' + self.__class__.taskNamesList[
    #             self.__class__.taskIndex] + '" data prepareing error: "'
    # + path + '" missed, no files deleted')
    #     self.testErrCodes.append(1)
    #     return

'''


""" ПЕРЕПИСАТЬ ФУНКИЦИЮ УНИВЕРСАЛЬНО """
def find_all_files(path_dir, pattern='', ext=''):
    """
    Returns a list with information about files.

    :param path_dir: path to a directory with any nesting
    :param pattern :
    :param ext     : file extension to find
    :return        : tuple(full_file_path, dir_file, file) where
                        full_file_path - absolute file path
                        dir_file - file directory
                        file     - file name
    """

    list_all_path_to_files = list()

    for dir_path, dir_name, list_filename in os.walk(path_dir):
        if ext:
            files = filter(
               lambda file_ext: file_ext if file_ext.endswith(ext) else False,
               list_filename
            )

        if pattern:
            files = filter(
               lambda file_ext: file_ext if pattern.lower() == file_ext.lower()
               else False, list_filename
            )

        for file in files:
            full_file_path = os.path.normpath(os.path.join(dir_path, file))

            if os.path.isfile(full_file_path):
                dir_file = os.path.split(os.path.dirname(full_file_path))[1]
                list_all_path_to_files.append((full_file_path, dir_file, file))

    return list_all_path_to_files

def parsing_file_variance(log_file):

    calculated_variance = {}
    with open(log_file, 'r') as file_var:
        list_str_from_remlf90 = list(map(
            lambda x: x.strip(), file_var.readlines()
        ))

        for item in list_str_from_remlf90:
            index = list_str_from_remlf90.index(item) + 1

            if 'Genetic variance(s)' in item:
                calculated_variance['genetic_variance'] = \
                    list_str_from_remlf90[index]

            elif 'Residual variance(s)' in item:
                calculated_variance['residual_variance'] = \
                    list_str_from_remlf90[index]

            elif '-2logL' and 'AIC' in item:
                div_item_on_part = item.split(':')

                for i_part in div_item_on_part:
                    if '-2logL' in i_part:
                        calculated_variance['BIC'] = \
                            i_part.split('=')[-1].strip()

                    elif 'AIC' in i_part:
                        calculated_variance['AIC'] = \
                            i_part.split('=')[-1].strip()

            else:
                pass

        calculated_variance['heritability'] = round(
            float(calculated_variance['genetic_variance']) / (
                    float(calculated_variance['genetic_variance']) +
                    float(calculated_variance['residual_variance'])
            ), 3)

    return calculated_variance



if __name__ == "__main__":

    path_to_file_model_name = './Input_model_data/list_model_name.json'
    path_to_file = './Input_model_data'

    with open(path_to_file_model_name, 'r') as file:
        model_name = json.loads(file.read())

    info_by_feature = {}
    for item_file in os.listdir(path_to_file):
        path_to_feature_dir = re_path(os.path.join(path_to_file, item_file))

        lst_info_by_model = []
        if os.path.isdir(path_to_feature_dir):
            lst_file = find_all_files(path_to_feature_dir,
                                     pattern='airemlf90.log')

            for item_data in lst_file:
                info_by_file = {}
                path_to_file_var, dir_name, file_name = item_data
                data_log = parsing_file_variance(path_to_file_var)

                info_by_file['model_name'] = model_name[dir_name]
                info_by_file.update(data_log)

                lst_info_by_model.append(info_by_file)


            info_by_feature[item_file] = lst_info_by_model


    with open("INFO_BY_MODEL.json", "w") as write_file:
        json.dump(info_by_feature, write_file, indent=4)


